import numpy as np
import networkx as nx
from scipy.sparse import coo_matrix, block_diag, coo_array
from .utils import SolutionTemporal

class TemporalSystemBuilder():
    def __init__(self,coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,inductive_mutuals_coords,inductive_mutuals_data,res_mutual_coords,res_mutual_data):
        """Class for building the linear system to solve for both temporal and frequency studies

        Parameters
        ----------
        coil_coords : np.array of shape (L,2) L being the number of coils
            Connectivity coordinates in the graph of each inductance
            Repetitions will be considered as another coil in //
        coil_data : np.array of shape (L,) L being the number of coils
            Values of inductance whose coordinates are in coil_coords
        res_coords : np.array of shape (R,2) R being the number of resistances
            Connectivity coordinates in the graph of each resistance
            Repetitions will be considered as another resistance in //
        res_data : np.array of shape (R,) R being the number of resistances
            Values of resistance whose coordinates are in res_coords
        capa_coords : np.array of shape (C,2) C being the number of capacities
            Connectivity coordinates in the graph of each capacity
            Repetitions will be considered as another capacity in //
        capa_data : np.array of shape (R,) R being the number of capacities
            Values of capacity whose coordinates are in res_coords
        inductive_mutuals_coords : np.array of shape (M,2) M being the number of mutuals in the system
            The coordinates are the indices of the coils in the coil_data array that have a mutual effect
            Inductive mutuals are only supported between coils
            Repetitions will be considered as another mutual
        inductive_mutuals_data : np.array of shape (M,) M being the number of mutuals in the system
            Values of mutuals whose coordinates are in inductive_mutuals_coords.
            The mutual follows the sign of the nodes given in the coil_coords array
        res_mutual_coords :np.array of shape (N,2) N being the number of resistive mutuals in the system
            the coordinates are the indices of the coils or resistance in the coil_data and res_data array that have a resistive mutual effect
            resitive mutuals are only supported between coils and resistances
        res_mutual_data : np.array of shape (N,) N being the number of resistive mutuals in the system
            Values of mutuals whose coordinates are in res_mutuals_coords.
            The mutual follows the sign of the nodes given in the coil_coords or res_coord array
        """

        self.coil_coords=coil_coords
        self.coil_data=coil_data
        self.res_coords=res_coords
        self.res_data=res_data
        self.capa_coords=capa_coords
        self.capa_data=capa_data
        self.inductive_mutuals_coords=inductive_mutuals_coords
        self.inductive_mutuals_data=inductive_mutuals_data
        self.res_mutual_coords=res_mutual_coords
        self.res_mutual_data=res_mutual_data
        self.current_sources_coords=np.zeros((2,0),dtype=int)
        self.current_sources_data=np.array([],dtype=int)
        self.voltage_sources_coords=np.zeros((2,0),dtype=int)
        self.voltage_sources_data=np.array([],dtype=int)
        self.source_count = 0

        ## initializing second member as empty
        self.rhs = (np.array([]),(np.array([],dtype=int),))

        self.analysed=False

    def graph_analysis(self):
        """Function to run the graph analysis:
        To be able to build the system one needs to figure out how many subsystems are present in the problem
        We test it with networkx and build some internal variables that are necessary for the system building
        In case you want to change the connectivity graph of your system you need to rerun the analysis to make check wether the topology changed or not
        """

        self.all_coords = np.concatenate((self.coil_coords,self.res_coords,self.capa_coords,self.voltage_sources_coords),axis=1)
        all_points = np.unique(self.all_coords)
        if all_points.shape != np.max(self.all_coords)+1:
            print("Warning: There is one or multiple lonely nodes please clean your impedence graph")
        self.all_impedences = np.concatenate([self.coil_data,self.res_data,self.capa_data,self.voltage_sources_data],axis=0)

        # actual number of node in the system
        self.size = np.max(self.all_coords)+1
        # number of intensities
        self.number_intensities = self.all_impedences.shape[0]
        ## making graph and checking number of subgraphs for ground enforcing and intensity checking
        unique_coords = np.unique(self.all_coords,axis=1)
        sym_graph = np.concatenate((unique_coords,np.stack((unique_coords[1],unique_coords[0]),axis=0)),axis=1)
        links = np.ones(sym_graph.shape[1])
        self.graph =  nx.from_scipy_sparse_array(coo_matrix((links,(sym_graph[0],sym_graph[1]))))
        ## keep the subgraphs
        self.list_of_subgraphs = [ list(sub) for sub in nx.connected_components(self.graph)]
        self.number_of_subsystems = len(self.list_of_subgraphs)
        ## location of ground
        self.affected_potentials = [-1]*self.number_of_subsystems
        ## by default remove 1 node equation per subsytem otherwise system is singular
        self.deleted_equation_current = [subsystem[0] for subsystem in self.list_of_subgraphs]
        ## shifter for intensities equations
        rescaler = np.zeros(self.size)
        rescaler[self.deleted_equation_current]=1
        rescaler = -np.cumsum(rescaler)
        self.rescaler =rescaler.astype(int)
        ## offsets for simplifying building the system
        offset_j = self.all_impedences.shape[0]
        offset_i = self.size-len(self.deleted_equation_current)
        self.offset_i = offset_i
        self.offset_j = offset_j

        ## State -> analysed
        self.analysed=True

    def set_ground(self,*args):
        """Function to affect a ground to subsystems
        If the system already has a ground provided then a warning is displayed and ground reaffected
        """
        ## Running graph analysis if not done
        if not self.analysed:
            self.graph_analysis()
        for index in args:
            for pivot,subsystem in enumerate(self.list_of_subgraphs):
                if index in subsystem:
                    if self.affected_potentials[pivot]!=-1:
                        print(f"Subsystem {pivot} already add a ground, reaffecting the value")
                    self.affected_potentials[pivot]=index
                    break

    def affect_potentials(self):
        """Function to check whether the grounds were all affected and assign some if some are missing
        """
        ## Running graph analysis if not done
        if not self.analysed:
            self.graph_analysis()
        for i in range(len(self.affected_potentials)):
            if -1 == self.affected_potentials[i]:
                self.affected_potentials[i]= self.list_of_subgraphs[i][0]
                print(f"Subsytem {i} has not been affected to the ground, we chose {self.list_of_subgraphs[i][0]}")

    def build_system(self):
        """Building sytem assuming that data was given as coords / data tuples
        This function builds 3 matrixes:
        S1 which is the real part of the system
        S2 which is the derivative part of the system
        S_init which is the system that needs to be solved for having initial conditions departing from null conditions
        """
        ## Running graph analysis if not done
        if not self.analysed:
            self.graph_analysis()
        ## affecting grounds if need be
        self.affect_potentials()
        ## Building rhs
        self.build_second_member()
        ## Building a sparse COO matrix


        ## Building all vectorized values necessary
        i_s_vals = np.max(self.all_coords,axis=0)
        j_s_vals = np.min(self.all_coords,axis=0)
        values = self.all_impedences


        ## node laws (only add to S1)
        data_nodes = np.concatenate((np.ones(self.number_intensities),-np.ones(self.number_intensities)),axis=0)
        j_s_nodes = np.tile(np.arange(self.number_intensities,dtype=int),(2,))
        i_s_nodes = np.concatenate((i_s_vals,j_s_vals),axis=0)

        # Removing one current equation per subsytem
        mask_removed_eq = ~np.isin(i_s_nodes,self.deleted_equation_current)
        data_nodes_S1 = data_nodes[mask_removed_eq]
        j_s_nodes_S1 = j_s_nodes[mask_removed_eq]
        i_s_nodes_S1 = i_s_nodes[mask_removed_eq]
        i_s_nodes_S1 = i_s_nodes_S1 + self.rescaler[i_s_nodes_S1]



        ## Kirchoff

        ## coils (contribution to S1 for potentials, S2 for intensities)
        i_s_coils = np.max(self.coil_coords,axis=0)
        j_s_coils = np.min(self.coil_coords,axis=0)
        i_s_edges_coil_S1 = self.offset_i + np.concatenate([np.arange(self.coil_data.shape[0],dtype=int)]*2,axis=0 )
        j_s_edges_coil_S1 = np.concatenate([self.offset_j+i_s_coils,self.offset_j+j_s_coils],axis=0)
        data_edges_coil_S1 = np.concatenate([np.ones(self.coil_data.shape[0]),-np.ones(self.coil_data.shape[0])],axis=0)

        i_s_edges_coil_S2 = self.offset_i + np.arange(self.coil_data.shape[0],dtype=int)
        j_s_edges_coil_S2 = np.arange(self.coil_data.shape[0],dtype=int)
        data_edges_coil_S2 = self.coil_data


        offset_coil= self.coil_data.shape[0]

        ## resistance (only goes to S1)
        i_s_res = np.max(self.res_coords,axis=0)
        j_s_res = np.min(self.res_coords,axis=0)
        i_s_edges_res_S1 = self.offset_i+offset_coil + np.concatenate([np.arange(self.res_data.shape[0],dtype=int)]*3,axis=0 )
        j_s_edges_res_S1 = np.concatenate([self.offset_j+i_s_res,self.offset_j+j_s_res,offset_coil+np.arange(self.res_data.shape[0],dtype=int)],axis=0)
        data_edges_res_S1 = np.concatenate([np.ones(self.res_data.shape[0]),-np.ones(self.res_data.shape[0]),self.res_data],axis=0)

        offset_res = self.res_data.shape[0]



        ## capacities (contribution to S2 for potentials, S1 for intensities)
        i_s_capa = np.max(self.capa_coords,axis=0)
        j_s_capa = np.min(self.capa_coords,axis=0)
        i_s_edges_capa_S1 = self.offset_i+offset_coil+offset_res + np.arange(self.capa_data.shape[0],dtype=int)
        j_s_edges_capa_S1 = offset_coil+offset_res+np.arange(self.capa_data.shape[0],dtype=int)
        data_edges_capa_S1 = np.ones(self.capa_data.shape[0])

        i_s_edges_capa_S2 = self.offset_i+offset_coil+offset_res + np.concatenate([np.arange(self.capa_data.shape[0],dtype=int)]*2,axis=0 )
        j_s_edges_capa_S2 = np.concatenate([self.offset_j+i_s_capa,self.offset_j+j_s_capa],axis=0)
        data_edges_capa_S2 = np.concatenate([self.capa_data,-self.capa_data],axis=0)


        ## adding mutuals to the system
        sign = np.sign(self.all_coords[0]-self.all_coords[1])

        ## inductive mutuals (contribution to S2 only)
        i_s_additionnal_S2 = self.offset_i + np.concatenate((self.inductive_mutuals_coords[0],self.inductive_mutuals_coords[1]),axis=0)
        j_s_additionnal_S2 = np.concatenate((self.inductive_mutuals_coords[1],self.inductive_mutuals_coords[0]),axis=0)
        data_additionnal_S2 = np.tile(self.inductive_mutuals_data*sign[self.inductive_mutuals_coords[0]]*sign[self.inductive_mutuals_coords[1]],(2,))

        ## resistive mutuals (contribution to S1 only)
        i_s_additionnal_S1= self.offset_i + np.concatenate((self.res_mutual_coords[0],self.res_mutual_coords[1]),axis=0)
        j_s_additionnal_S1 = np.concatenate((self.res_mutual_coords[1],self.res_mutual_coords[0]),axis=0)
        data_additionnal_S1 = np.tile(self.res_mutual_data*sign[self.res_mutual_coords[0]]*sign[self.res_mutual_coords[1]],(2,))


        ## ground equations (1 per subsystem)
        i_s_ground_S1 = np.arange(self.offset_i+self.offset_j,self.size+self.offset_j) - self.source_count
        j_s_ground_S1 = self.offset_j+np.array(self.affected_potentials)
        data_ground_S1 = np.ones(len(self.affected_potentials))

        ## voltage sources equations
        i_s_sources_S1 = np.tile(np.arange(0,self.voltage_sources_data.shape[0]),2)+self.number_intensities+self.size - self.source_count
        j_s_sources_S1 = np.concatenate((self.offset_j + self.voltage_sources_coords[0],self.offset_j+self.voltage_sources_coords[1]),axis=0)
        data_source_S1 = np.concatenate((np.ones_like(self.voltage_sources_data),-np.ones_like(self.voltage_sources_data)))


        ## building S1 system
        i_s_S1 = np.concatenate((i_s_nodes_S1,i_s_edges_coil_S1,i_s_edges_res_S1,i_s_edges_capa_S1,i_s_additionnal_S1,i_s_ground_S1,i_s_sources_S1),axis=0)
        j_s_S1 = np.concatenate((j_s_nodes_S1,j_s_edges_coil_S1,j_s_edges_res_S1,j_s_edges_capa_S1,j_s_additionnal_S1,j_s_ground_S1,j_s_sources_S1),axis=0)
        data_S1 = np.concatenate((data_nodes_S1,data_edges_coil_S1,data_edges_res_S1,data_edges_capa_S1,data_additionnal_S1,data_ground_S1,data_source_S1),axis=0)

        ## building S2 system
        i_s_S2 = np.concatenate((i_s_edges_coil_S2,i_s_edges_capa_S2,i_s_additionnal_S2),axis=0)
        j_s_S2 = np.concatenate((j_s_edges_coil_S2,j_s_edges_capa_S2,j_s_additionnal_S2),axis=0)
        data_S2 = np.concatenate((data_edges_coil_S2,data_edges_capa_S2,data_additionnal_S2),axis=0)

        ## building S_init system
        i_s_init = np.concatenate((i_s_nodes_S1,i_s_edges_coil_S2,i_s_edges_res_S1,i_s_edges_capa_S2,i_s_ground_S1,i_s_sources_S1),axis=0)
        j_s_init = np.concatenate((j_s_nodes_S1,j_s_edges_coil_S2,j_s_edges_res_S1,j_s_edges_capa_S2,j_s_ground_S1,j_s_sources_S1),axis=0)
        data_init = np.concatenate((data_nodes_S1,data_edges_coil_S2,data_edges_res_S1,data_edges_capa_S2,data_ground_S1,data_source_S1),axis=0)

        self.S_init=(data_init,(i_s_init,j_s_init))

        self.S1 = (data_S1,(i_s_S1.astype(int),j_s_S1.astype(int)))
        self.S2 = (data_S2,(i_s_S2.astype(int),j_s_S2.astype(int)))
        ## return for debug if needed
        return self.S1,self.S2,self.S_init

    def build_second_member(self,check=True):
        """General second member builder, need to be called after graph analysis and voltage and current sources tests

        Parameters
        ----------
        check : bool, optional
            Whether to check or not that current injections append on the same subsystem, setting to False speeds up, by default True

        Returns
        -------
        tuple
            rhs tuple in case it is needed

        Raises
        ------
        IndexError
            raised when current injection does not belong to the same subsystem
        """
        ## Testing current
        if check:
            for i in range(self.current_sources_coords.shape[1]):
                input_node=self.current_sources_coords[0,i]
                output_node=self.current_sources_coords[1,i]
                if self.number_of_subsystems>=2:
                    valid = False
                    for system in self.list_of_subgraphs:
                        if input_node in system and output_node in system:
                            valid =True
                        else:
                            continue
                    if not valid:
                        raise IndexError("Nodes indicated do not belong to the same subsystem")


        ## Building current injection
        in_current_nodes = self.current_sources_coords[0]
        in_current_data = - self.current_sources_data
        out_current_nodes = self.current_sources_coords[1]
        out_current_data = self.current_sources_data
        current_nodes = np.concatenate((in_current_nodes,out_current_nodes),axis=0)
        current_data = np.concatenate((in_current_data,out_current_data),axis=0)

        # Removing deleted node equations
        mask_removed_eq = ~np.isin(current_nodes,self.deleted_equation_current)
        current_nodes = current_nodes[mask_removed_eq]
        current_data = current_data[mask_removed_eq]
        current_nodes = current_nodes + self.rescaler[current_nodes]

        ## Building voltage rhs
        voltage_nodes = np.arange(0,self.voltage_sources_data.shape[0])+self.number_intensities+self.size - self.source_count
        voltage_data = self.voltage_sources_data

        self.rhs=(np.concatenate((current_data,voltage_data),axis=0),(np.concatenate((current_nodes,voltage_nodes),axis=0),))
        return self.rhs

    def add_current_source(self,intensity,input_node,output_node):
        """Function to build a second member for the scenario of current injection in the sparse system

        Parameters
        ----------
        intensity : float
            intensity to inject
        input_node : int
            which node to take for injection
        output_node : int
            which node for current retrieval
        """
        self.current_sources_coords = np.append(self.current_sources_coords,np.array([[input_node],[output_node]]),axis=1)
        self.current_sources_data = np.append(self.current_sources_data,np.array([intensity]))

    def add_voltage_source(self,voltage,input_node,output_node):
        """Adding a voltage source to the system. This adds one equation and one degree of freedom in the system (the source intensity)

        Parameters
        ----------
        voltage : complex
            enforced voltage
        input_node : int
            node where the voltage is enforced
        output_node : int
            node from where the voltage is enforced

        """
        self.voltage_sources_coords = np.append(self.voltage_sources_coords,np.array([[input_node],[output_node]]),axis=1)
        self.voltage_sources_data = np.append(self.voltage_sources_data,np.array([voltage]))
        self.source_count+=1

    def get_init_system(self,sparse_rhs=False):
        """Function to determine the initial state of the system as coo_matrix
        Parameters
        ----------
        sparse_rhs : bool
            whether to return the second member as a sparse vector

        Returns
        -------
        sys: coo_matrix
            left hand side of the equation to solve
        rhs:
            right hand side of the system to solve
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

    def get_system(self,sparse_rhs=False):
        """Function to get the whole system to solve as coo_matrix
        Parameters
        ----------
        sparse_rhs : bool
            whether to return the second member as a sparse vector

        Returns
        -------
        sys1 : coo_matrix
            real part of the system
        sys2 : coo_matrix
            img part of the system (need to be multiplyied by j.omega for frequency studies)
        rhs: np.array
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

    def get_frequency_system(self,omega,sparse_rhs=False):
        """Function to get the complex matrix for frequency studies

        Parameters
        ----------
        omega : float
            pulsation of the system
        sparse_rhs : bool
            whether to return the second member as a sparse vector

        Returns
        -------
        sys: coo_matrix, dtype=complex
            complex system for frequency studies
        rhs: np.array
            right hand side of the system
        """
        sys1,sys2,rhs = self.get_system(sparse_rhs=sparse_rhs)
        return sys1+1j*omega*sys2,rhs


    def build_intensity_and_voltage_from_vector(self,sol):
        """Utility function to reshape the solution of frequency studies or temporal studies

        Parameters
        ----------
        sol : np.array of shape (*, self.number_intensities+self.size+self.source_count)
            array containing as many solutions as wanted

        Returns
        -------
        coil_intensities: np.array of shape (*,coil_data.shape[0])
            intensities following the directions given by coil_coords
        res_intensities: np.array of shape (*,res_data.shape[0])
            intensities following the directions given by res_coords
        capa_intensities: np.array of shape (*,capa_data.shape[0])
            intensities following the directions given by res_coords
        voltages: np.array of shape (*,self.size)
            voltage of the points in the graph
        source_intensities: np.array of shape (*, self.source_count)
            if there are voltage sources this contains the source intensities
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

