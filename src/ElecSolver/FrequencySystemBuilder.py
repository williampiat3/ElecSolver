import numpy as np
import networkx as nx
from scipy.sparse import coo_matrix


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
        all_points = np.unique(impedence_coords)
        if all_points.shape != np.max(impedence_coords)+1:
            IndexError("There is one or multiple lonely nodes please clean your impedence graph")
        self.impedence_coords = impedence_coords
        self.impedence_data = impedence_data
        self.mutuals_coords = mutuals_coords
        self.mutuals_data = mutuals_data
        self.size = np.max(impedence_coords)+1
        self.number_intensities = self.impedence_data.shape[0]
        ## making graph and checking number of subgraphs
        unique_coords = np.unique(self.impedence_coords,axis=1)
        sym_graph = np.concatenate((unique_coords,np.stack((unique_coords[1],unique_coords[0]),axis=0)),axis=1)
        links = np.ones(sym_graph.shape[1])
        self.graph =  nx.from_scipy_sparse_array(coo_matrix((links,(sym_graph[0],sym_graph[1]))))
        self.list_of_subgraphs = [ list(sub) for sub in nx.connected_components(self.graph)]
        self.number_of_subsystems = len(self.list_of_subgraphs)
        self.affected_potentials = [-1]*self.number_of_subsystems
        self.deleted_equation_current = [subsystem[0] for subsystem in self.list_of_subgraphs]
        ## number of tension sources
        self.source_count = 0
        self.source_signs = np.array([],dtype=int)
        self.rhs = (np.array([]),(np.array([],dtype=int),))
        rescaler = np.zeros(self.size)
        rescaler[self.deleted_equation_current]=1
        rescaler = -np.cumsum(rescaler)
        self.rescaler =rescaler.astype(int)
        offset_j = self.impedence_data.shape[0]
        offset_i = self.size-len(self.deleted_equation_current)
        self.offset_i = offset_i
        self.offset_j = offset_j

    def set_mass(self,*args):
        for index in args:
            for pivot,subsystem in enumerate(self.list_of_subgraphs):
                if index in subsystem:
                    if self.affected_potentials[pivot]!=-1:
                        print(f"Subsystem {pivot} already add a mass, reaffecting the value")
                    self.affected_potentials[pivot]=index
                    break

    def affect_potentials(self):
        """Function to check whether the masses were affected and assign some if not
        """
        for i in range(len(self.affected_potentials)):
            if -1 == self.affected_potentials[i]:
                self.affected_potentials[i]= self.list_of_subgraphs[i][0]
                print(f"Subsytem {i} has not been affected to the mass, we chose {self.list_of_subgraphs[i][0]}")

    def build_system(self):
        """Building sytem assuming that data was given as COO matrices
        Fully vectorized for best performance
        it is faster but less understandable
        """
        ## affecting masses if need be
        self.affect_potentials()
        ## Building a sparse COO matrix


        ## Building all vectorized values necessary
        i_s_vals = np.max(self.impedence_coords,axis=0)
        j_s_vals = np.min(self.impedence_coords,axis=0)
        values = self.impedence_data


        ## node laws
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
        i_s_edges = self.offset_i + np.concatenate([np.arange(self.number_intensities,dtype=int)]*3,axis=0 )
        j_s_edges = np.concatenate([self.offset_j+i_s_vals,self.offset_j+j_s_vals,np.arange(self.number_intensities,dtype=int)],axis=0)
        data_edges = np.concatenate([np.ones(values.shape[0],dtype=complex),-np.ones(values.shape[0],dtype=complex),values],axis=0)



        ## adding mutuals to the system
        sign = np.sign(self.impedence_coords[0]-self.impedence_coords[1])


        i_s_additionnal = self.offset_i + np.concatenate((self.mutuals_coords[0],self.mutuals_coords[1]),axis=0)
        j_s_additionnal = np.concatenate((self.mutuals_coords[1],self.mutuals_coords[0]),axis=0)
        data_additionnal = np.tile(self.mutuals_data*sign[self.mutuals_coords[0]]*sign[self.mutuals_coords[1]],(2,))


        ## mass equations (1 per subsystem)
        i_s_mass = np.arange(self.offset_i+self.offset_j,self.offset_i+self.offset_j+len(self.affected_potentials))
        j_s_mass = self.offset_j+np.array(self.affected_potentials)
        data_mass = np.ones(len(self.affected_potentials))


        i_s = np.concatenate((i_s_nodes,i_s_edges,i_s_additionnal,i_s_mass),axis=0)
        j_s = np.concatenate((j_s_nodes,j_s_edges,j_s_additionnal,j_s_mass),axis=0)
        data = np.concatenate((data_nodes,data_edges,data_additionnal,data_mass),axis=0)

        self.system = (data,(i_s.astype(int),j_s.astype(int)))
        return self.system

    def get_system(self):
        """Function to get the system
        Returns
        -------
        sys: scipy.coo_matrix
            Linear system to solve
        rhs: np.ndarray
            Second member of the system
        """
        (data_rhs,(nodes,)) = self.rhs
        (data,(i_s,j_s)) = self.system
        sys = coo_matrix(self.system)
        rhs = np.zeros(self.number_intensities+self.size+self.source_count)
        np.add.at(rhs, nodes, data_rhs)
        return sys,rhs

    def build_second_member_intensity(self,intensity,input_node,output_node):
        """Function to build a second member for the scenario of current injection in the sparse system

        Parameters
        ----------
        intensity : float
            intensity to inject
        input_node : int
            which node to take for injection
        output_node : int
            which node for current retrieval

        Returns
        -------
        Tuple(np.array,(np.array))
            second member of the linear system in a COO like format
        """
        data = []
        nodes = []
        if self.number_of_subsystems>=2:
            valid = False
            for system in self.list_of_subgraphs:
                if input_node in system and output_node in system:
                    valid =True
                else:
                    continue
            if not valid:
                raise IndexError("Nodes indicated do not belong to the same subsystem")

        if input_node not in self.deleted_equation_current:
            nodes.append(input_node+self.rescaler[input_node])
            data.append(-intensity) # inbound intensity
        if output_node not in self.deleted_equation_current:
            nodes.append(output_node+self.rescaler[output_node])
            data.append(intensity) # outbound intensity
        (data_rhs,(nodes_rhs,)) = self.rhs
        self.rhs = (np.append(data_rhs,np.array(data),axis=0),(np.append(nodes_rhs,np.array(nodes),axis=0),))
        return self.rhs

    def build_second_member_tension(self,tension,input_node,output_node):
        """building second member in the case of an enforced tension
        To be able to enforce a tension we need to remove one of the equation in the system
        and add the tension enforced equation
        By default we remove the first kirchoff law within the correct subsystem
        Warning this function changes the system, please re-get the system before solving it (call get_system)

        Parameters
        ----------
        tension : complex
            enforced tension
        input_node : int
            node where the tension is enforce
        output_node : int
            node from where the tension is enforced

        Returns
        -------
        Tuple(np.array,(np.array))
            second member of the linear system in a COO like format
        """
        targeted_subsystem = 0
        if self.number_of_subsystems>=2:
            valid = False
            for index_sub,system in enumerate(self.list_of_subgraphs):
                if input_node in system and output_node in system:
                    valid =True
                    targeted_subsystem=index_sub
                else:
                    continue
            if not valid:
                raise IndexError("Nodes indicated do not belong to the same subsystem")
        print("Warning this function changes the system, please re-evaluate it before solving the system (call get_system)")
        (data,(i_s,j_s)) = self.system
        (data_rhs,(nodes,)) = self.rhs

        new_eq = self.number_intensities+self.size +self.source_count
        sign = np.sign(output_node-input_node)

        ## Adding the new equation in S1
        if input_node in self.deleted_equation_current:
            data = np.append(data,np.array([1.,-1.,-sign]),axis=0)
            i_s = np.append(i_s,np.array([new_eq,new_eq,output_node+self.rescaler[output_node]]),axis=0)
            j_s = np.append(j_s,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq]),axis=0)
        elif output_node in self.deleted_equation_current:
            data = np.append(data,np.array([1.,-1.,sign]),axis=0)
            i_s = np.append(i_s,np.array([new_eq,new_eq,input_node+self.rescaler[input_node]]),axis=0)
            j_s = np.append(j_s,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq]),axis=0)
        else:
            data = np.append(data,np.array([1.,-1.,sign,-sign]),axis=0)
            i_s = np.append(i_s,np.array([new_eq,new_eq,input_node+self.rescaler[input_node],output_node+self.rescaler[output_node]]),axis=0)
            j_s = np.append(j_s,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq,new_eq]),axis=0)

        ## second member value
        nodes= np.append(nodes, [new_eq],axis=0)
        data_rhs= np.append(data_rhs, [tension],axis=0)

         ## reaffecting the systems and second member
        self.system = (data,(i_s,j_s))
        self.rhs = (data_rhs,(nodes,))

        ## counter for tension source
        self.source_count +=1
        self.source_signs=np.append(self.source_signs,[sign])
        return self.rhs

    def build_intensity_and_voltage_from_vector(self,sol):
        sign = np.sign(self.impedence_coords[1]-self.impedence_coords[0])
        if self.source_count!=0:
            return (sol[:self.number_intensities]*sign,
                    sol[self.number_intensities:-self.source_count],
                    sol[...,-self.source_count:]*self.source_signs
                    )

        else:
            return (sol[:self.number_intensities]*sign,
                    sol[self.number_intensities:],
                    np.array([],dtype=float)
                    )