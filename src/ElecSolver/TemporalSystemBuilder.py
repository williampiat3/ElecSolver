import numpy as np
from scipy.sparse import coo_matrix, block_diag, coo_array
from .utils import SolutionTemporal,GradientsParametersTemporal,compute_graph_components
import warnings

class TemporalSystemBuilder():
    def __init__(self,coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,inductive_mutual_coords,inductive_mutual_data,res_mutual_coords,res_mutual_data):
        """Initialize a temporal electrical system builder.

        Parameters
        ----------
        coil_coords : numpy.ndarray, shape (2, L)
            Node coordinates for each inductor. Repeated coordinates represent
            parallel inductors.
        coil_data : numpy.ndarray, shape (L,)
            Inductance values corresponding to ``coil_coords``.
        res_coords : numpy.ndarray, shape (2, R)
            Node coordinates for each resistor. Repeated coordinates represent
            parallel resistors.
        res_data : numpy.ndarray, shape (R,)
            Resistance values corresponding to ``res_coords``.
        capa_coords : numpy.ndarray, shape (2, C)
            Node coordinates for each capacitor. Repeated coordinates represent
            parallel capacitors.
        capa_data : numpy.ndarray, shape (C,)
            Capacitance values corresponding to ``capa_coords``.
        inductive_mutual_coords : numpy.ndarray, shape (2, M)
            Inductor-index pairs with mutual inductive coupling.
        inductive_mutual_data : numpy.ndarray, shape (M,)
            Mutual inductance values. Each sign follows the node order in
            ``coil_coords``.
        res_mutual_coords : numpy.ndarray, shape (2, P)
            Component-index pairs with mutual resistive coupling. Resistive
            mutuals are supported between inductors and resistors.
        res_mutual_data : numpy.ndarray, shape (P,)
            Mutual resistance values. Each sign follows the node order in
            ``coil_coords`` or ``res_coords``.
        """

        self.coil_coords=coil_coords
        self.coil_data=coil_data
        self.res_coords=res_coords
        self.res_data=res_data
        self.capa_coords=capa_coords
        self.capa_data=capa_data
        self.inductive_mutual_coords=inductive_mutual_coords
        self.inductive_mutual_data=inductive_mutual_data
        self.res_mutual_coords=res_mutual_coords
        self.res_mutual_data=res_mutual_data
        # Source arrays are populated after initialization so they can be validated.
        self.current_source_coords=np.zeros((2,0),dtype=int)
        self.current_source_data=np.array([],dtype=int)
        self.voltage_source_coords=np.zeros((2,0),dtype=int)
        self.voltage_source_data=np.array([],dtype=int)
        self.source_count = 0

        # Initialize the right-hand side as empty.
        self.rhs = (np.array([]),(np.array([],dtype=int),))

        self.analysed=False

    def graph_analysis(self):
        """Analyze connected components and initialize system dimensions.

        Rerun this method after changing the connectivity graph so that topology
        changes are reflected in the assembled system.

        Raises
        ------
        RuntimeError
            If no voltage or current source has been added.
        """
        if self.voltage_source_data.shape==0 and self.current_source_data.shape==0:
            raise RuntimeError("The system does not have any sources (voltage or current) defined before running the graph analysis. Please define all your sources before calling set_ground, build_system or graph_analysis")

        self.all_coords = np.concatenate((self.coil_coords,self.res_coords,self.capa_coords,self.voltage_source_coords),axis=1)
        all_points = np.unique(self.all_coords)
        if all_points.shape != np.max(self.all_coords)+1:
            warnings.warn("Warning: There is one or multiple lonely nodes please clean your impedence graph")

        if self.analysed:
            warnings.warn("Warning: analysis was already performed: grounds will be reasigned")

        self.all_impedences = np.concatenate([self.coil_data,self.res_data,self.capa_data,self.voltage_source_data],axis=0)

        # Actual number of nodes in the system.
        self.size = np.max(self.all_coords)+1
        # Number of currents.
        self.number_intensities = self.all_impedences.shape[0]
        # Keep the connected subgraphs.
        self.list_of_subgraphs = [ list(sub) for sub in compute_graph_components(self.all_coords)]
        self.number_of_subsystems = len(self.list_of_subgraphs)
        # Reassign previous grounds when repeating the analysis.
        if self.analysed:
            grounds_placeholder = self.affected_potentials
            self.affected_potentials = self.affected_potentials[:min(self.number_of_subsystems,len(self.affected_potentials))]
            self.set_ground(*grounds_placeholder)
        else:
            self.affected_potentials = [-1]*self.number_of_subsystems

        # Remove one node equation per subsystem to avoid a singular system.
        self.deleted_equation_current = [subsystem[0] for subsystem in self.list_of_subgraphs]
        # Compute offsets for current-equation indices to account for the removed node equations.
        rescaler = np.zeros(self.size)
        rescaler[self.deleted_equation_current]=1
        rescaler = -np.cumsum(rescaler)
        self.rescaler =rescaler.astype(int)
        # Offsets used while building the system.
        offset_j = self.all_impedences.shape[0]
        offset_i = self.size-len(self.deleted_equation_current)
        self.offset_i = offset_i
        self.offset_j = offset_j

        # Mark the builder as analyzed.
        self.analysed=True

    def get_nx_graph(self):
        """Return the NetworkX graph representation of the system.

        NetworkX is imported locally to avoid making it a module-level dependency.

        Returns
        -------
        networkx.Graph
            Graph representation of the system.
        """
        import networkx as nx
        # Run analysis if needed so ``all_coords`` is defined once.
        if not self.analysed:
            self.graph_analysis()

        unique_coords = np.unique(self.all_coords,axis=1)
        sym_graph = np.concatenate((unique_coords,np.stack((unique_coords[1],unique_coords[0]),axis=0)),axis=1)
        links = np.ones(sym_graph.shape[1])
        graph =  nx.from_scipy_sparse_array(coo_matrix((links,(sym_graph[0],sym_graph[1]))))
        return graph


    def set_ground(self,*args):
        """Assign ground nodes to subsystems.

        Parameters
        ----------
        *args : int
            Node indices to assign as grounds.
        """
        # Run graph analysis if needed.
        if not self.analysed:
            self.graph_analysis()
        for index in args:
            for pivot,subsystem in enumerate(self.list_of_subgraphs):
                if index in subsystem:
                    if self.affected_potentials[pivot]!=-1 and self.affected_potentials[pivot]!=index:
                        print(f"Subsystem {pivot} already add a ground, reaffecting the value")
                    self.affected_potentials[pivot]=index
                    break

    def affect_potentials(self):
        """Assign default ground nodes to subsystems without explicit grounds."""
        # Run graph analysis if needed.
        if not self.analysed:
            self.graph_analysis()
        for i in range(len(self.affected_potentials)):
            if -1 == self.affected_potentials[i]:
                self.affected_potentials[i]= self.list_of_subgraphs[i][0]
                print(f"Subsystem {i} has not been assigned a ground; using {self.list_of_subgraphs[i][0]}")

    def build_system(self):
        """Build the real, derivative, and initial-condition systems.

        This function builds three matrices:
        - S1 which is the real part of the system (invertible)
        - S2 which is the derivative part of the system (singular)
        - S_init which is the system that needs to be solved for having initial conditions (invertible)

        Let ``N`` be the number of nodes, ``M`` the number of impedances, ``k``
        the number of subsystems, and ``s`` the number of voltage sources. The
        unknown vector contains ``M + s`` currents followed by ``N`` potentials.
        Let rhs be the second member of the system containing current injections and voltage sources contributions.
        The temporal system is then given by:
        U_0 = S_init^{-1}.rhs
        S1.U + S2.\\dot{U} = rhs

        The equations are ordered as follows:

        - node laws (N-k equations)
        - impedance equations (M equations)
        - ground equations (k equations)
        - voltage sources equations (s equations)

        Returns
        -------
        S1 : tuple
            COO data and coordinates for the real system.
        S2 : tuple
            COO data and coordinates for the derivative system.
        S_init : tuple
            COO data and coordinates for the initial-condition system.
        """
        # Run graph analysis if needed.
        if not self.analysed:
            self.graph_analysis()
        # Assign grounds if needed.
        self.affect_potentials()
        # Build the right-hand side.
        self.build_second_member()

        # Start building the system

        # Build the vectorized COO entries.
        i_s_vals = np.max(self.all_coords,axis=0)
        j_s_vals = np.min(self.all_coords,axis=0)
        values = self.all_impedences


        # Node laws contribute only to S1.
        data_nodes = np.concatenate((np.ones(self.number_intensities),-np.ones(self.number_intensities)),axis=0)
        j_s_nodes = np.tile(np.arange(self.number_intensities,dtype=int),(2,))
        i_s_nodes = np.concatenate((i_s_vals,j_s_vals),axis=0)

        # Remove one current equation per subsystem.
        mask_removed_eq = ~np.isin(i_s_nodes,self.deleted_equation_current)
        data_nodes_S1 = data_nodes[mask_removed_eq]
        j_s_nodes_S1 = j_s_nodes[mask_removed_eq]
        i_s_nodes_S1 = i_s_nodes[mask_removed_eq]
        i_s_nodes_S1 = i_s_nodes_S1 + self.rescaler[i_s_nodes_S1]


        # Impedance equations.

        # Inductors contribute potentials to S1 and currents to S2.
        i_s_coils = np.max(self.coil_coords,axis=0)
        j_s_coils = np.min(self.coil_coords,axis=0)
        i_s_edges_coil_S1 = self.offset_i + np.concatenate([np.arange(self.coil_data.shape[0],dtype=int)]*2,axis=0 )
        j_s_edges_coil_S1 = np.concatenate([self.offset_j+i_s_coils,self.offset_j+j_s_coils],axis=0)
        data_edges_coil_S1 = np.concatenate([np.ones(self.coil_data.shape[0]),-np.ones(self.coil_data.shape[0])],axis=0)

        i_s_edges_coil_S2 = self.offset_i + np.arange(self.coil_data.shape[0],dtype=int)
        j_s_edges_coil_S2 = np.arange(self.coil_data.shape[0],dtype=int)
        data_edges_coil_S2 = self.coil_data


        offset_coil= self.coil_data.shape[0]

        # Resistors contribute only to S1.
        i_s_res = np.max(self.res_coords,axis=0)
        j_s_res = np.min(self.res_coords,axis=0)
        i_s_edges_res_S1 = self.offset_i+offset_coil + np.concatenate([np.arange(self.res_data.shape[0],dtype=int)]*3,axis=0 )
        j_s_edges_res_S1 = np.concatenate([self.offset_j+i_s_res,self.offset_j+j_s_res,offset_coil+np.arange(self.res_data.shape[0],dtype=int)],axis=0)
        data_edges_res_S1 = np.concatenate([np.ones(self.res_data.shape[0]),-np.ones(self.res_data.shape[0]),self.res_data],axis=0)

        offset_res = self.res_data.shape[0]



        # Capacitors contribute potentials to S2 and currents to S1.
        i_s_capa = np.max(self.capa_coords,axis=0)
        j_s_capa = np.min(self.capa_coords,axis=0)
        i_s_edges_capa_S1 = self.offset_i+offset_coil+offset_res + np.arange(self.capa_data.shape[0],dtype=int)
        j_s_edges_capa_S1 = offset_coil+offset_res+np.arange(self.capa_data.shape[0],dtype=int)
        data_edges_capa_S1 = np.ones(self.capa_data.shape[0])

        i_s_edges_capa_S2 = self.offset_i+offset_coil+offset_res + np.concatenate([np.arange(self.capa_data.shape[0],dtype=int)]*2,axis=0 )
        j_s_edges_capa_S2 = np.concatenate([self.offset_j+i_s_capa,self.offset_j+j_s_capa],axis=0)
        data_edges_capa_S2 = np.concatenate([self.capa_data,-self.capa_data],axis=0)


        # Add mutual couplings to the system.
        sign = np.sign(self.all_coords[0]-self.all_coords[1])

        # Inductive mutuals contribute only to S2.
        i_s_additionnal_S2 = self.offset_i + np.concatenate((self.inductive_mutual_coords[0],self.inductive_mutual_coords[1]),axis=0)
        j_s_additionnal_S2 = np.concatenate((self.inductive_mutual_coords[1],self.inductive_mutual_coords[0]),axis=0)
        data_additionnal_S2 = np.tile(self.inductive_mutual_data*sign[self.inductive_mutual_coords[0]]*sign[self.inductive_mutual_coords[1]],(2,))

        # Resistive mutuals contribute only to S1.
        i_s_additionnal_S1= self.offset_i + np.concatenate((self.res_mutual_coords[0],self.res_mutual_coords[1]),axis=0)
        j_s_additionnal_S1 = np.concatenate((self.res_mutual_coords[1],self.res_mutual_coords[0]),axis=0)
        data_additionnal_S1 = np.tile(self.res_mutual_data*sign[self.res_mutual_coords[0]]*sign[self.res_mutual_coords[1]],(2,))


        # Add one ground equation per subsystem.
        i_s_ground_S1 = np.arange(self.offset_i+self.offset_j,self.size+self.offset_j) - self.source_count
        j_s_ground_S1 = self.offset_j+np.array(self.affected_potentials)
        data_ground_S1 = np.ones(len(self.affected_potentials))

        # Add voltage-source equations.
        i_s_sources_S1 = np.tile(np.arange(0,self.voltage_source_data.shape[0]),2)+self.number_intensities+self.size - self.source_count
        j_s_sources_S1 = np.concatenate((self.offset_j + self.voltage_source_coords[0],self.offset_j+self.voltage_source_coords[1]),axis=0)
        data_source_S1 = np.concatenate((np.ones_like(self.voltage_source_data),-np.ones_like(self.voltage_source_data)))


        # Build the S1 system.
        i_s_S1 = np.concatenate((i_s_nodes_S1,i_s_edges_coil_S1,i_s_edges_res_S1,i_s_edges_capa_S1,i_s_additionnal_S1,i_s_ground_S1,i_s_sources_S1),axis=0)
        j_s_S1 = np.concatenate((j_s_nodes_S1,j_s_edges_coil_S1,j_s_edges_res_S1,j_s_edges_capa_S1,j_s_additionnal_S1,j_s_ground_S1,j_s_sources_S1),axis=0)
        data_S1 = np.concatenate((data_nodes_S1,data_edges_coil_S1,data_edges_res_S1,data_edges_capa_S1,data_additionnal_S1,data_ground_S1,data_source_S1),axis=0)

        # Build reverse S1 indices for gradients.
        self.S1_reverse_res = data_nodes_S1.shape[0]+data_edges_coil_S1.shape[0] + self.res_data.shape[0]*2 + np.arange(self.res_data.shape[0])
        self.S1_reverse_res_sign = np.ones_like(self.res_data)
        # resistive mutuals are used twice in the system but only contribute once to the gradient
        self.S1_reverse_mutual_res_1 = data_nodes_S1.shape[0]+data_edges_coil_S1.shape[0] + data_edges_res_S1.shape[0] + data_edges_capa_S1.shape[0] + np.arange(self.res_mutual_data.shape[0])
        self.S1_reverse_mutual_res_2 = data_nodes_S1.shape[0]+data_edges_coil_S1.shape[0] + data_edges_res_S1.shape[0] + data_edges_capa_S1.shape[0] + np.arange(self.res_mutual_data.shape[0],self.res_mutual_data.shape[0]*2)
        self.S1_reverse_mutual_res_sign_1 = sign[self.res_mutual_coords[0]]*sign[self.res_mutual_coords[1]]
        self.S1_reverse_mutual_res_sign_2 = self.S1_reverse_mutual_res_sign_1

        # Build the S2 system.
        i_s_S2 = np.concatenate((i_s_edges_coil_S2,i_s_edges_capa_S2,i_s_additionnal_S2),axis=0)
        j_s_S2 = np.concatenate((j_s_edges_coil_S2,j_s_edges_capa_S2,j_s_additionnal_S2),axis=0)
        data_S2 = np.concatenate((data_edges_coil_S2,data_edges_capa_S2,data_additionnal_S2),axis=0)

        # Build reverse S2 indices for gradients.
        self.S2_reverse_coil = np.arange(data_edges_coil_S2.shape[0])
        self.S2_reverse_coil_sign = np.ones_like(data_edges_coil_S2)
        # Capacitors are used twice in the system but only contribute once to the gradient.
        self.S2_reverse_capa_1 = data_edges_coil_S2.shape[0] + np.arange(self.capa_data.shape[0])
        self.S2_reverse_capa_2 = data_edges_coil_S2.shape[0] + np.arange(self.capa_data.shape[0], self.capa_data.shape[0]*2)
        self.S2_reverse_capa_sign_1 = np.ones_like(self.capa_data)
        self.S2_reverse_capa_sign_2 = -np.ones_like(self.capa_data)
        # inductive mutuals are used twice in the system but only contribute once to the gradient
        self.S2_reverse_mutual_inductive_1 = data_edges_coil_S2.shape[0] + data_edges_capa_S2.shape[0] + np.arange(self.inductive_mutual_data.shape[0])
        self.S2_reverse_mutual_inductive_2 = data_edges_coil_S2.shape[0] + data_edges_capa_S2.shape[0] + np.arange(self.inductive_mutual_data.shape[0], self.inductive_mutual_data.shape[0]*2)
        self.S2_reverse_mutual_inductive_sign_1 = sign[self.inductive_mutual_coords[0]]*sign[self.inductive_mutual_coords[1]]
        self.S2_reverse_mutual_inductive_sign_2 = self.S2_reverse_mutual_inductive_sign_1

        # Build the initial-condition system.
        i_s_init = np.concatenate((i_s_nodes_S1,i_s_edges_coil_S2,i_s_edges_res_S1,i_s_edges_capa_S2,i_s_ground_S1,i_s_sources_S1),axis=0)
        j_s_init = np.concatenate((j_s_nodes_S1,j_s_edges_coil_S2,j_s_edges_res_S1,j_s_edges_capa_S2,j_s_ground_S1,j_s_sources_S1),axis=0)
        data_init = np.concatenate((data_nodes_S1,data_edges_coil_S2,data_edges_res_S1,data_edges_capa_S2,data_ground_S1,data_source_S1),axis=0)

        # Only resistors contribute to initial-condition system gradients.
        self.S_init_reverse_res = data_nodes_S1.shape[0]+data_edges_coil_S2.shape[0] + self.res_data.shape[0]*2 + np.arange(self.res_data.shape[0])
        self.S_init_reverse_res_sign = np.ones_like(self.res_data)

        self.S_init=(data_init,(i_s_init,j_s_init))

        self.S1 = (data_S1,(i_s_S1.astype(int),j_s_S1.astype(int)))
        self.S2 = (data_S2,(i_s_S2.astype(int),j_s_S2.astype(int)))
        # Return tuple forms for debugging or direct use.
        return self.S1,self.S2,self.S_init

    def build_second_member(self,check=True):
        """Build the right-hand side after graph and source validation.

        Parameters
        ----------
        check : bool, optional
            Whether to verify that each current injection connects nodes in the
            same subsystem.

        Returns
        -------
        tuple
            COO data and coordinates for the right-hand side.

        Raises
        ------
        IndexError
            If a current injection connects different subsystems.
        """
        # Validate current sources.
        if check:
            for i in range(self.current_source_coords.shape[1]):
                input_node=self.current_source_coords[0,i]
                output_node=self.current_source_coords[1,i]
                if self.number_of_subsystems>=2:
                    valid = False
                    for system in self.list_of_subgraphs:
                        if input_node in system and output_node in system:
                            valid =True
                        else:
                            continue
                    if not valid:
                        raise IndexError(f"Nodes {input_node} and {output_node} do not belong to the same subsystem, can't create a current source between these two points")


        # Build current injections.
        in_current_nodes = self.current_source_coords[0]
        in_current_data = - self.current_source_data
        out_current_nodes = self.current_source_coords[1]
        out_current_data = self.current_source_data
        current_nodes = np.concatenate((in_current_nodes,out_current_nodes),axis=0)
        current_data = np.concatenate((in_current_data,out_current_data),axis=0)

        # Removing deleted node equations
        mask_removed_eq = ~np.isin(current_nodes,self.deleted_equation_current)
        current_nodes = current_nodes[mask_removed_eq]
        current_data = current_data[mask_removed_eq]
        current_nodes = current_nodes + self.rescaler[current_nodes]

        # Build voltage-source terms.
        voltage_nodes = np.arange(0,self.voltage_source_data.shape[0])+self.number_intensities+self.size - self.source_count
        voltage_data = self.voltage_source_data

        self.rhs=(np.concatenate((current_data,voltage_data),axis=0),(np.concatenate((current_nodes,voltage_nodes),axis=0),))
        return self.rhs

    def add_current_source(self,intensity,input_node,output_node):
        """Add a current source to the system.

        Parameters
        ----------
        intensity : float
            Current to inject.
        input_node : int
            Injection node.
        output_node : int
            Current retrieval node.
        """
        self.current_source_coords = np.append(self.current_source_coords,np.array([[input_node],[output_node]]),axis=1)
        self.current_source_data = np.append(self.current_source_data,np.array([intensity]))

    def add_voltage_source(self,voltage,input_node,output_node):
        """Add a voltage source and its current degree of freedom.

        Parameters
        ----------
        voltage : complex
            Enforced voltage.
        input_node : int
            Node where the voltage is enforced.
        output_node : int
            Reference node for the enforced voltage.
        """
        if self.analysed:
            warnings.warn("Warning: adding a tension source when analysis is performed may result in system topology change. You may need to rerun graph_analysis if it is the case.")
        self.voltage_source_coords = np.append(self.voltage_source_coords,np.array([[input_node],[output_node]]),axis=1)
        self.voltage_source_data = np.append(self.voltage_source_data,np.array([voltage]))
        self.source_count+=1

    def get_init_system(self,sparse_rhs=True):
        """Return the initial-condition system and right-hand side.

        Parameters
        ----------
        sparse_rhs : bool, optional
            Whether to return the right-hand side as a sparse COO array.

        Returns
        -------
        sys : scipy.sparse.coo_matrix, shape (n, n)
            Initial-condition system matrix.
        rhs : scipy.sparse.coo_array or numpy.ndarray, shape (n,)
            Right-hand side of the system.
        """
        size = self.number_intensities+self.size
        sys = coo_matrix(self.S_init,shape=(size,size))
        if sparse_rhs:
            (data_rhs,(nodes,)) = self.rhs
            rhs = coo_array((data_rhs,(nodes,)),shape=(size,))
            rhs.sum_duplicates()
        else:
            rhs = np.zeros(size)
            (data_rhs,(nodes,)) = self.rhs
            np.add.at(rhs, nodes, data_rhs)
        return sys,rhs

    def get_system(self,sparse_rhs=True):
        """Return the real and derivative systems and right-hand side.

        Parameters
        ----------
        sparse_rhs : bool, optional
            Whether to return the right-hand side as a sparse COO array.

        Returns
        -------
        sys1 : scipy.sparse.coo_matrix, shape (n, n)
            Real part of the system.
        sys2 : scipy.sparse.coo_matrix, shape (n, n)
            Derivative part of the system. For frequency studies, multiply it
            by ``1j * omega``.
        rhs : scipy.sparse.coo_array or numpy.ndarray, shape (n,)
            Right-hand side of the system.
        """
        size = self.number_intensities+self.size
        sys1 = coo_matrix(self.S1,shape=(size,size))
        sys2 = coo_matrix(self.S2,shape=(size,size))
        if sparse_rhs:
            (data_rhs,(nodes,)) = self.rhs
            rhs = coo_array((data_rhs,(nodes,)),shape=(size,))
            rhs.sum_duplicates()
        else:
            rhs = np.zeros(size)
            (data_rhs,(nodes,)) = self.rhs
            np.add.at(rhs, nodes, data_rhs)
        return sys1,sys2,rhs

    def backpropagate_gradients(self, dS1=None, dS2=None, dS_init=None, drhs=None):
        """Backpropagate system gradients to electrical parameters.

        Call this method after :meth:`build_system`. Gradients correspond to
        sparse data arrays, not full dense matrices or vectors. Contributions
        from repeated parameter values are summed.

        Parameters
        ----------
        dS1 : numpy.ndarray, optional
            Gradient of the S1 data array.
        dS2 : numpy.ndarray, optional
            Gradient of the S2 data array.
        dS_init : numpy.ndarray, optional
            Gradient of the initial-condition system data array.
        drhs : numpy.ndarray, optional
            Gradient of the right-hand-side data array.

        Returns
        -------
        GradientsParametersTemporal
            Gradients with shapes matching their corresponding parameter arrays.
        """
        # Initialize parameter gradients.
        grads_coils = np.zeros_like(self.coil_data,dtype=float)
        grads_res = np.zeros_like(self.res_data,dtype=float)
        grads_capa = np.zeros_like(self.capa_data,dtype=float)
        grads_voltage_sources = np.zeros_like(self.voltage_source_data,dtype=float)
        grads_inductive_mutuals = np.zeros_like(self.inductive_mutual_data,dtype=float)
        grads_res_mutuals = np.zeros_like(self.res_mutual_data,dtype=float)
        grads_current_sources = np.zeros_like(self.current_source_data,dtype=float)

        if dS1 is not None:
            # Resistor contribution to S1.
            grads_res += dS1[self.S1_reverse_res]*self.S1_reverse_res_sign
            # Resistive mutual contributions to S1.
            grads_res_mutuals += dS1[self.S1_reverse_mutual_res_1]*self.S1_reverse_mutual_res_sign_1
            grads_res_mutuals += dS1[self.S1_reverse_mutual_res_2]*self.S1_reverse_mutual_res_sign_2


        if dS2 is not None:
            # Inductor contribution to S2.
            grads_coils += dS2[self.S2_reverse_coil]*self.S2_reverse_coil_sign
            # Capacitor contribution to S2.
            grads_capa += dS2[self.S2_reverse_capa_1]*self.S2_reverse_capa_sign_1
            grads_capa += dS2[self.S2_reverse_capa_2]*self.S2_reverse_capa_sign_2
            # Inductive mutual contributions to S2.
            grads_inductive_mutuals += dS2[self.S2_reverse_mutual_inductive_1]*self.S2_reverse_mutual_inductive_sign_1
            grads_inductive_mutuals += dS2[self.S2_reverse_mutual_inductive_2]*self.S2_reverse_mutual_inductive_sign_2

        if dS_init is not None:
            grads_res += dS_init[self.S_init_reverse_res]*self.S_init_reverse_res_sign

        if drhs is not None:
            # Handle current sources separately because some equations are removed.
            in_current_nodes = self.current_source_coords[0]
            out_current_nodes = self.current_source_coords[1]
            in_current_sign = -np.ones_like(in_current_nodes)
            out_current_sign = np.ones_like(out_current_nodes)
            mask_eq_in = ~np.isin(in_current_nodes,self.deleted_equation_current)
            mask_eq_out = ~np.isin(out_current_nodes,self.deleted_equation_current)
            in_current_nodes = in_current_nodes[mask_eq_in]
            in_current_nodes = in_current_nodes + self.rescaler[in_current_nodes]
            in_current_sign = in_current_sign[mask_eq_in]
            out_current_nodes = out_current_nodes[mask_eq_out]
            out_current_nodes = out_current_nodes + self.rescaler[out_current_nodes]
            out_current_sign = out_current_sign[mask_eq_out]

            grads_current_sources[mask_eq_in] += drhs[np.arange(in_current_nodes.shape[0])]*in_current_sign
            grads_current_sources[mask_eq_out] += drhs[in_current_nodes.shape[0] + np.arange(out_current_nodes.shape[0])]*out_current_sign
            grads_voltage_sources += drhs[(in_current_nodes.shape[0] + out_current_nodes.shape[0]):]
        return GradientsParametersTemporal(grads_coils, grads_res, grads_capa, grads_inductive_mutuals, grads_res_mutuals, grads_voltage_sources, grads_current_sources)


    def get_frequency_system(self,omega,sparse_rhs=False):
        """Return the complex system for a given ``omega`` (angular frequency).

        Parameters
        ----------
        omega : float
            Angular frequency of the system.
        sparse_rhs : bool, optional
            Whether to return the right-hand side as a sparse COO array.

        Returns
        -------
        sys : scipy.sparse.coo_matrix, shape (n, n)
            Complex system for frequency studies.
        rhs : scipy.sparse.coo_array or numpy.ndarray, shape (n,)
            Right-hand side of the system.
        """
        sys1,sys2,rhs = self.get_system(sparse_rhs=sparse_rhs)
        return sys1+1j*omega*sys2,rhs


    def build_intensity_and_voltage_from_vector(self,sol):
        """Convert a raw solution vector to named solution components.

        Parameters
        ----------
        sol : numpy.ndarray, shape (..., self.number_intensities + self.size)
            One or more raw solution vectors.

        Returns
        -------
        SolutionTemporal
            Inductor, resistor, and capacitor currents; node potentials; and
            voltage-source currents.
        """
        sign = np.sign(self.all_coords[1]-self.all_coords[0])
        offset_coil = self.coil_data.shape[0]
        offset_res = self.res_data.shape[0]
        offset_capa = self.capa_data.shape[0]
        if self.source_count!=0:
            return SolutionTemporal(sol[...,:offset_coil]*sign[:offset_coil],
                sol[...,offset_coil:offset_coil+offset_res]*sign[offset_coil:offset_coil+offset_res],
                sol[...,offset_coil+offset_res:offset_coil+offset_res+offset_capa]*sign[offset_coil+offset_res:offset_coil+offset_res+offset_capa],
                sol[...,self.number_intensities:],
                sol[...,offset_coil+offset_res+offset_capa:offset_coil+offset_res+offset_capa+self.source_count]*sign[offset_coil+offset_res+offset_capa:offset_coil+offset_res+offset_capa+self.source_count]

            )
        else:
            return SolutionTemporal(sol[...,:offset_coil]*sign[:offset_coil],
                sol[...,offset_coil:offset_coil+offset_res]*sign[offset_coil:offset_coil+offset_res],
                sol[...,offset_coil+offset_res:offset_coil+offset_res+offset_capa]*sign[offset_coil+offset_res:offset_coil+offset_res+offset_capa],
                sol[...,self.number_intensities:],
                np.array([],dtype=float)
            )

    def build_vector_from_intensity_and_voltage(self,solution:SolutionTemporal):
        """Build a raw solution vector from named solution components.

        Parameters
        ----------
        solution : SolutionTemporal
            Components of the solution.

        Returns
        -------
        numpy.ndarray, shape (..., self.number_intensities + self.size)
            One or more raw solution vectors.
        """
        sign = np.sign(self.all_coords[1]-self.all_coords[0])
        offset_coil = self.coil_data.shape[0]
        offset_res = self.res_data.shape[0]
        offset_capa = self.capa_data.shape[0]
        if self.source_count!=0:
            return np.concatenate((solution.coil_intensities*sign[:offset_coil],
                solution.res_intensities*sign[offset_coil:offset_coil+offset_res],
                solution.capa_intensities*sign[offset_coil+offset_res:offset_coil+offset_res+offset_capa],
                solution.voltages,
                solution.source_intensities*sign[offset_coil+offset_res+offset_capa:offset_coil+offset_res+offset_capa+self.source_count]
            ),axis=solution.voltages.ndim-1)
        else:
            return np.concatenate((solution.coil_intensities*sign[:offset_coil],
                solution.res_intensities*sign[offset_coil:offset_coil+offset_res],
                solution.capa_intensities*sign[offset_coil+offset_res:offset_coil+offset_res+offset_capa],
                solution.voltages
            ),axis=solution.voltages.ndim-1)

