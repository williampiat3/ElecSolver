import numpy as np
from scipy.sparse import coo_matrix,coo_array
from .utils import SolutionFrequency,GradientsParametersFrequency,compute_graph_components
import warnings


class FrequencySystemBuilder():
    def __init__(self,impedence_coords,impedence_data,mutual_coords,mutual_data):
        """Initialize a frequency-domain electrical system builder.

        The builder supports complex values and general impedances.

        Parameters
        ----------
        impedence_coords : numpy.ndarray, shape (2, N)
            Node coordinates for each impedance.
        impedence_data : numpy.ndarray, shape (N,)
            Impedance values between nodes ``impedence_coords[:, i]``.
        mutual_coords : numpy.ndarray, shape (2, M)
            Indices of impedance pairs with mutual coupling.
        mutual_data : numpy.ndarray, shape (M,)
            Mutual values between the indexed impedance pairs. Each sign follows
            the node order in ``impedence_coords``.
        """
        self.impedence_coords = impedence_coords
        self.impedence_data = impedence_data
        self.mutual_coords = mutual_coords
        self.mutual_data = mutual_data
        self.current_source_coords=np.zeros((2,0),dtype=int)
        self.current_source_data=np.array([],dtype=int)
        self.voltage_source_coords=np.zeros((2,0),dtype=int)
        self.voltage_source_data=np.array([],dtype=int)
        self.source_count = 0
        # Initialize the right-hand side as empty.
        self.rhs = (np.array([]),(np.array([],dtype=int),))

        self.analysed=False

    def graph_analysis(self):
        """Analyze connected components and initialize system dimensions."""
        self.all_coords = np.concatenate((self.impedence_coords,self.voltage_source_coords),axis=1)
        all_points = np.unique(self.all_coords)
        if all_points.shape != np.max(self.all_coords)+1:
            warnings.warn("Warning: There is one or multiple lonely nodes please clean your impedence graph")

        if self.analysed:
            warnings.warn("Warning: analysis was already performed: grounds will be reasigned")

        self.all_impedences = np.concatenate([self.impedence_data,self.voltage_source_data],axis=0)

        # Actual number of nodes in the system.
        self.size = np.max(self.all_coords)+1
        # Number of currents.
        self.number_intensities = self.all_impedences.shape[0]
        # Keep the connected subgraphs.
        self.list_of_subgraphs = [ sub.tolist() for sub in compute_graph_components(self.all_coords)]
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
        # rescaler is an offset for current-equation indices after removing equations.
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
        """Build the complex system in COO format.

        Let ``N`` be the number of nodes, ``M`` the number of impedances, ``k``
        the number of subsystems, and ``s`` the number of voltage sources. The
        system has size ``N + M + s``. Its unknowns are the ``M + s`` currents,
        followed by the ``N`` potentials.

        The equations are ordered as follows:

        - node laws (N-k equations)
        - Kirchhoff laws (M equations)
        - ground equations (k equations)
        - voltage sources equations (s equations)

        Returns
        -------
        tuple
            COO data and coordinate arrays for the system matrix.
        """
        # Run graph analysis if needed.
        if not self.analysed:
            self.graph_analysis()
        # Assign grounds if needed.
        self.affect_potentials()
        # Build the right-hand side.
        self.build_second_member()

        # Starting to build the system matrix.

        # Build node-law entries.
        i_s_vals = np.max(self.all_coords,axis=0)
        j_s_vals = np.min(self.all_coords,axis=0)
        data_nodes = np.concatenate((np.ones(self.number_intensities),-np.ones(self.number_intensities)),axis=0)
        j_s_nodes = np.tile(np.arange(self.number_intensities,dtype=int),(2,))
        i_s_nodes =np.concatenate((i_s_vals,j_s_vals),axis=0)


        # Remove one current equation per subsystem.
        mask_removed_eq = ~np.isin(i_s_nodes,self.deleted_equation_current)
        data_nodes = data_nodes[mask_removed_eq]
        j_s_nodes = j_s_nodes[mask_removed_eq]
        i_s_nodes = i_s_nodes[mask_removed_eq]
        i_s_nodes = i_s_nodes + self.rescaler[i_s_nodes]



        # Build impedence entries.
        i_s_vals = np.max(self.impedence_coords,axis=0)
        j_s_vals = np.min(self.impedence_coords,axis=0)
        values = self.impedence_data
        i_s_edges = self.offset_i + np.concatenate([np.arange(self.impedence_data.shape[0],dtype=int)]*3,axis=0 )
        j_s_edges = np.concatenate([self.offset_j+i_s_vals,self.offset_j+j_s_vals,np.arange(self.impedence_data.shape[0],dtype=int)],axis=0)
        data_edges = np.concatenate([np.ones(values.shape[0],dtype=complex),-np.ones(values.shape[0],dtype=complex),values],axis=0)



        # Add mutual couplings to the system.
        sign = np.sign(self.impedence_coords[0]-self.impedence_coords[1])


        i_s_additionnal = self.offset_i + np.concatenate((self.mutual_coords[0],self.mutual_coords[1]),axis=0)
        j_s_additionnal = np.concatenate((self.mutual_coords[1],self.mutual_coords[0]),axis=0)
        data_additionnal = np.tile(self.mutual_data*sign[self.mutual_coords[0]]*sign[self.mutual_coords[1]],(2,))


        # Add one ground equation per subsystem.
        i_s_ground = np.arange(self.offset_i+self.offset_j,self.size+self.offset_j) - self.source_count
        j_s_ground = self.offset_j+np.array(self.affected_potentials)
        data_ground = np.ones(len(self.affected_potentials))

        # Add voltage-source equations.
        i_s_sources = np.tile(np.arange(0,self.voltage_source_data.shape[0]),2)+self.number_intensities+self.size - self.source_count
        j_s_sources = np.concatenate((self.offset_j + self.voltage_source_coords[0],self.offset_j+self.voltage_source_coords[1]),axis=0)
        data_source = np.concatenate((np.ones_like(self.voltage_source_data),-np.ones_like(self.voltage_source_data)))

        # Build the sparse system.
        i_s = np.concatenate((i_s_nodes,i_s_edges,i_s_additionnal,i_s_ground,i_s_sources),axis=0)
        j_s = np.concatenate((j_s_nodes,j_s_edges,j_s_additionnal,j_s_ground,j_s_sources),axis=0)
        data = np.concatenate((data_nodes,data_edges,data_additionnal,data_ground,data_source),axis=0)

        # Build reverse system indices for gradients.
        self.S_reverse_impedence = data_nodes.shape[0]+self.impedence_data.shape[0]*2 + np.arange(self.impedence_data.shape[0])
        self.S_reverse_impedence_sign = np.ones_like(self.impedence_data)
        # mutuals are used twice in the system but only contribute once to the gradient
        self.S_reverse_mutual_1 = data_nodes.shape[0]+data_edges.shape[0] + np.arange(self.mutual_data.shape[0])
        self.S_reverse_mutual_2 = data_nodes.shape[0]+data_edges.shape[0] + np.arange(self.mutual_data.shape[0],2*self.mutual_data.shape[0])
        self.S_reverse_mutual_sign_1 = sign[self.mutual_coords[0]]*sign[self.mutual_coords[1]]
        self.S_reverse_mutual_sign_2 = sign[self.mutual_coords[0]]*sign[self.mutual_coords[1]]

        self.system = (data,(i_s.astype(int),j_s.astype(int)))
        return self.system

    def build_second_member(self,check=True):
        """Build the right-hand side after graph and source validation.

        Parameters
        ----------
        check : bool, optional
            Whether to verify that each current injection connects nodes in the
            same subsystem. Setting to False will speed up a little bit.

        Returns
        -------
        tuple
            COO data and coordinates for the right-hand side in the form of a tuple (data, (nodes,)).

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

    def get_system(self,sparse_rhs=True):
        """Return the assembled system matrix and right-hand side.

        Parameters
        ----------
        sparse_rhs : bool, optional
            Whether to return the right-hand side as a sparse COO array.

        Returns
        -------
        sys : scipy.sparse.coo_matrix, shape (n, n)
            Linear system to solve.
        rhs : scipy.sparse.coo_array or numpy.ndarray, shape (n,)
            Right-hand side of the system.
        """
        size = self.number_intensities+self.size
        (data_rhs,(nodes,)) = self.rhs
        sys = coo_matrix(self.system,shape=(size,size))
        if sparse_rhs:
            (data_rhs,(nodes,)) = self.rhs
            rhs = coo_array((data_rhs,(nodes,)),shape=(size,))
            rhs.sum_duplicates()
        else:
            rhs = np.zeros(size)
            (data_rhs,(nodes,)) = self.rhs
            np.add.at(rhs, nodes, data_rhs)
        return sys,rhs

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
        if self.analysed == True:
            warnings.warn("Warning: adding a tension source when analysis is performed may result in system topology change. You may need to rerun graph_analysis if it is the case.")
        self.voltage_source_coords = np.append(self.voltage_source_coords,np.array([[input_node],[output_node]]),axis=1)
        self.voltage_source_data = np.append(self.voltage_source_data,np.array([voltage]))
        self.source_count+=1

    def backpropagate_gradients(self, dS=None, drhs=None):
        """Backpropagate system gradients to electrical parameters.

        Call this method after :meth:`build_system`. Gradients correspond to
        the sparse data arrays, not full dense matrices or vectors. Contributions
        from repeated parameter values are summed.

        Parameters
        ----------
        dS : numpy.ndarray, optional
            Gradient of the system matrix data array.
        drhs : numpy.ndarray, optional
            Gradient of the right-hand-side data array.

        Returns
        -------
        GradientsParametersFrequency
            Gradients with shapes matching their corresponding parameter arrays.
        """
        # Initialize parameter gradients.
        grads_impedence = np.zeros_like(self.impedence_data,dtype=complex)
        grads_voltage_sources = np.zeros_like(self.voltage_source_data,dtype=complex)
        grads_mutual = np.zeros_like(self.mutual_data,dtype=complex)
        grads_current_sources = np.zeros_like(self.current_source_data,dtype=complex)

        if dS is not None:
            # Impedance contribution to the system.
            grads_impedence += dS[self.S_reverse_impedence]*self.S_reverse_impedence_sign
            # Mutual contributions to the system.
            grads_mutual += dS[self.S_reverse_mutual_1]*self.S_reverse_mutual_sign_1
            grads_mutual += dS[self.S_reverse_mutual_2]*self.S_reverse_mutual_sign_2


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
        return GradientsParametersFrequency(grads_impedence, grads_mutual, grads_voltage_sources, grads_current_sources)



    def build_intensity_and_voltage_from_vector(self,sol):
        """Convert a raw solution vector to named solution components.

        Parameters
        ----------
        sol : numpy.ndarray, shape (..., n)
            Raw system solution.

        Returns
        -------
        SolutionFrequency
            Impedance currents, node potentials, and voltage-source currents.
        """
        sign = np.sign(self.all_coords[1]-self.all_coords[0])
        if self.source_count!=0:
            return SolutionFrequency(sol[...,:self.number_intensities-self.source_count]*sign[:self.number_intensities-self.source_count],
                    sol[...,self.number_intensities:],
                    sol[...,self.number_intensities-self.source_count:self.number_intensities]*sign[self.number_intensities-self.source_count:self.number_intensities]
                    )

        else:
            return SolutionFrequency(sol[...,:self.number_intensities]*sign,
                    sol[...,self.number_intensities:],
                    np.array([],dtype=float)
                    )
    def build_vector_from_intensity_and_voltage(self,solution: SolutionFrequency):
        """Build a raw solution vector from named solution components.

        Parameters
        ----------
        solution : SolutionFrequency
            Components of the solution.

        Returns
        -------
        numpy.ndarray, shape (..., self.number_intensities + self.size)
            One or more raw solution vectors.
        """
        sign = np.sign(self.all_coords[1]-self.all_coords[0])
        if self.source_count!=0:
            return np.concatenate((solution.intensities*sign[:self.number_intensities-self.source_count],
                    solution.potentials,
                    solution.intensities_sources*sign[self.number_intensities-self.source_count:self.number_intensities]
                    ),axis=solution.intensities.ndim-1)
        else:
            return np.concatenate((solution.intensities*sign,
                    solution.potentials
                    ),axis=solution.intensities.ndim-1)