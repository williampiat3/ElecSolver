import numpy as np
import networkx as nx
from scipy.sparse import coo_matrix,coo_array
from .utils import SolutionFrequency
import warnings


class FrequencySystemBuilder():
    def __init__(self,impedence_coords,impedence_data,mutuals_coords,mutuals_data):
        """FrequencySystemBuilder class for building an electrical sparse system
        that can be solved by any sparse solver
        it supports all forms of complex making this class fit for non linear impedences

        Parameters
        ----------
        impedence_coords : np.array of ints, shape = (2,N)
            impedence coordinates
        impedence_data : np.array of complex, shape = (N,)
            impedence value between points impedence_coords[:,i]
        mutuals_coords : np.array of ints, shape = (2,M)
            indexes of the impedence within impedence data which have a mutual
        mutuals_data : np.array of complex, shape=(M,)
            mutual value between impedences impedence_data[mutuals_coords[0,i]] impedence_data[mutuals_coords[1,i]]
            The mutual follows the order given in impedence coords
        """
        self.impedence_coords = impedence_coords
        self.impedence_data = impedence_data
        self.mutuals_coords = mutuals_coords
        self.mutuals_data = mutuals_data
        self.current_sources_coords=np.zeros((2,0),dtype=int)
        self.current_sources_data=np.array([],dtype=int)
        self.voltage_sources_coords=np.zeros((2,0),dtype=int)
        self.voltage_sources_data=np.array([],dtype=int)
        self.source_count = 0
        ## initializing second member as empty
        self.rhs = (np.array([]),(np.array([],dtype=int),))

        self.analysed=False

    def graph_analysis(self):
        self.all_coords = np.concatenate((self.impedence_coords,self.voltage_sources_coords),axis=1)
        all_points = np.unique(self.all_coords)
        if all_points.shape != np.max(self.all_coords)+1:
            warnings.warn("Warning: There is one or multiple lonely nodes please clean your impedence graph")

        if self.analysed:
            warnings.warn("Warning: analysis was already performed: grounds will be reasigned")

        self.all_impedences = np.concatenate([self.impedence_data,self.voltage_sources_data],axis=0)

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
        ## If analysis was already performed we take the previous grounds and try to reassign them to the system
        if self.analysed:
            grounds_placeholder = self.affected_potentials
            self.affected_potentials = self.affected_potentials[:min(self.number_of_subsystems,len(self.affected_potentials))]
            self.set_ground(*grounds_placeholder)
        else:
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
                    if self.affected_potentials[pivot]!=-1 and self.affected_potentials[pivot]!=index:
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
        """Building sytem assuming that data was given as COO matrices
        Fully vectorized for best performance
        it is faster but less understandable
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



        ## node laws
        i_s_vals = np.max(self.all_coords,axis=0)
        j_s_vals = np.min(self.all_coords,axis=0)
        data_nodes = np.concatenate((np.ones(self.number_intensities),-np.ones(self.number_intensities)),axis=0)
        j_s_nodes = np.tile(np.arange(self.number_intensities,dtype=int),(2,))
        i_s_nodes =np.concatenate((i_s_vals,j_s_vals),axis=0)


        # Removing one current equation per subsytem
        mask_removed_eq = ~np.isin(i_s_nodes,self.deleted_equation_current)
        data_nodes = data_nodes[mask_removed_eq]
        j_s_nodes = j_s_nodes[mask_removed_eq]
        i_s_nodes = i_s_nodes[mask_removed_eq]
        i_s_nodes = i_s_nodes + self.rescaler[i_s_nodes]



        ## Kirchoff
        i_s_vals = np.max(self.impedence_coords,axis=0)
        j_s_vals = np.min(self.impedence_coords,axis=0)
        values = self.impedence_data
        i_s_edges = self.offset_i + np.concatenate([np.arange(self.impedence_data.shape[0],dtype=int)]*3,axis=0 )
        j_s_edges = np.concatenate([self.offset_j+i_s_vals,self.offset_j+j_s_vals,np.arange(self.impedence_data.shape[0],dtype=int)],axis=0)
        data_edges = np.concatenate([np.ones(values.shape[0],dtype=complex),-np.ones(values.shape[0],dtype=complex),values],axis=0)



        ## adding mutuals to the system
        sign = np.sign(self.impedence_coords[0]-self.impedence_coords[1])


        i_s_additionnal = self.offset_i + np.concatenate((self.mutuals_coords[0],self.mutuals_coords[1]),axis=0)
        j_s_additionnal = np.concatenate((self.mutuals_coords[1],self.mutuals_coords[0]),axis=0)
        data_additionnal = np.tile(self.mutuals_data*sign[self.mutuals_coords[0]]*sign[self.mutuals_coords[1]],(2,))


        ## ground equations (1 per subsystem)
        i_s_ground = np.arange(self.offset_i+self.offset_j,self.size+self.offset_j) - self.source_count
        j_s_ground = self.offset_j+np.array(self.affected_potentials)
        data_ground = np.ones(len(self.affected_potentials))

        ## voltage sources equations
        i_s_sources = np.tile(np.arange(0,self.voltage_sources_data.shape[0]),2)+self.number_intensities+self.size - self.source_count
        j_s_sources = np.concatenate((self.offset_j + self.voltage_sources_coords[0],self.offset_j+self.voltage_sources_coords[1]),axis=0)
        data_source = np.concatenate((np.ones_like(self.voltage_sources_data),-np.ones_like(self.voltage_sources_data)))


        i_s = np.concatenate((i_s_nodes,i_s_edges,i_s_additionnal,i_s_ground,i_s_sources),axis=0)
        j_s = np.concatenate((j_s_nodes,j_s_edges,j_s_additionnal,j_s_ground,j_s_sources),axis=0)
        data = np.concatenate((data_nodes,data_edges,data_additionnal,data_ground,data_source),axis=0)

        self.system = (data,(i_s.astype(int),j_s.astype(int)))
        return self.system

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

    def get_system(self,sparse_rhs=False):
        """Function to get the system
        Parameters:
        -------
        sparse_rhs: bool, optionnal
            Whether to return the second member as a sparse.coo_array

        Returns
        -------
        sys: scipy.coo_matrix
            Linear system to solve
        rhs: np.ndarray
            Second member of the system
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
        if self.analysed == True:
            warnings.warn("Warning: adding a tension source when analysis is performed may result in system topology change. You may need to rerun graph_analysis if it is the case.")
        self.voltage_sources_coords = np.append(self.voltage_sources_coords,np.array([[input_node],[output_node]]),axis=1)
        self.voltage_sources_data = np.append(self.voltage_sources_data,np.array([voltage]))
        self.source_count+=1


    def build_intensity_and_voltage_from_vector(self,sol):
        sign = np.sign(self.all_coords[1]-self.all_coords[0])
        if self.source_count!=0:
            return SolutionFrequency(sol[:self.number_intensities-self.source_count]*sign[:self.number_intensities-self.source_count],
                    sol[self.number_intensities:],
                    sol[...,self.number_intensities-self.source_count:self.number_intensities]*sign[self.number_intensities-self.source_count:self.number_intensities]
                    )

        else:
            return SolutionFrequency(sol[:self.number_intensities]*sign,
                    sol[self.number_intensities:],
                    np.array([],dtype=float)
                    )