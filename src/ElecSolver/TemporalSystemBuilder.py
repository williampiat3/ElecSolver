import numpy as np
import networkx as nx
import sparse
from scipy.sparse import coo_matrix


class TemporalElectricSystemBuilder():
    def __init__(self,coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,inductive_mutuals_coords,inductive_mutuals_data,res_mutual_coords,res_mutual_data):
        """
        """
        self.all_coords = np.concatenate((coil_coords,res_coords,capa_coords),axis=1)
        all_points = np.unique(self.all_coords)
        if all_points.shape != np.max(self.all_coords)+1:
            IndexError("There is one or multiple lonely nodes please clean your impedence graph")
        self.all_impedences = np.concatenate([coil_data,res_data,capa_data],axis=0)

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



        self.size = np.max(self.all_coords)+1
        self.number_intensities = self.all_impedences.shape[0]
        ## making graph and checking number of subgraphs
        unique_coords = np.unique(self.all_coords,axis=1)
        sym_graph = np.concatenate((unique_coords,np.stack((unique_coords[1],unique_coords[0]),axis=0)),axis=1)
        links = np.ones(sym_graph.shape[1])
        self.graph =  nx.from_scipy_sparse_array(coo_matrix((links,(sym_graph[0],sym_graph[1]))))
        self.list_of_subgraphs = [ list(sub) for sub in nx.connected_components(self.graph)]
        self.number_of_subsystems = len(self.list_of_subgraphs)
        self.affected_potentials = [-1]*self.number_of_subsystems
        self.deleted_equation_current = [subsystem[0] for subsystem in self.list_of_subgraphs]
        self.source_count = 0
        self.source_signs = np.array([],dtype=int)
        self.rhs = (np.array([]),(np.array([],dtype=int),))
        rescaler = np.zeros(self.size)
        rescaler[self.deleted_equation_current]=1
        rescaler = -np.cumsum(rescaler)
        self.rescaler =rescaler.astype(int)
        offset_j = self.all_impedences.shape[0]
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
        i_s_vals = np.max(self.all_coords,axis=0)
        j_s_vals = np.min(self.all_coords,axis=0)
        values = self.all_impedences


        ## node laws (only to S1)
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

        ## coils
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



        ## capacities
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

        ## inductive mutuals
        i_s_additionnal_S2 = self.offset_i + np.concatenate((self.inductive_mutuals_coords[0],self.inductive_mutuals_coords[1]),axis=0)
        j_s_additionnal_S2 = np.concatenate((self.inductive_mutuals_coords[1],self.inductive_mutuals_coords[0]),axis=0)
        data_additionnal_S2 = np.tile(self.inductive_mutuals_data*sign[self.inductive_mutuals_coords[0]]*sign[self.inductive_mutuals_coords[1]],(2,))


        i_s_additionnal_S1= self.offset_i + np.concatenate((self.res_mutual_coords[0],self.res_mutual_coords[1]),axis=0)
        j_s_additionnal_S1 = np.concatenate((self.res_mutual_coords[1],self.res_mutual_coords[0]),axis=0)
        data_additionnal_S1 = np.tile(self.res_mutual_data*sign[self.res_mutual_coords[0]]*sign[self.res_mutual_coords[1]],(2,))


        ## mass equations (1 per subsystem)
        i_s_mass_S1 = np.arange(self.offset_i+self.offset_j,self.offset_i+self.offset_j+len(self.affected_potentials))
        j_s_mass_S1 = self.offset_j+np.array(self.affected_potentials)
        data_mass_S1 = np.ones(len(self.affected_potentials))




        i_s_S1 = np.concatenate((i_s_nodes_S1,i_s_edges_coil_S1,i_s_edges_res_S1,i_s_edges_capa_S1,i_s_additionnal_S1,i_s_mass_S1),axis=0)
        j_s_S1 = np.concatenate((j_s_nodes_S1,j_s_edges_coil_S1,j_s_edges_res_S1,j_s_edges_capa_S1,j_s_additionnal_S1,j_s_mass_S1),axis=0)
        data_S1 = np.concatenate((data_nodes_S1,data_edges_coil_S1,data_edges_res_S1,data_edges_capa_S1,data_additionnal_S1,data_mass_S1),axis=0)

        i_s_S2 = np.concatenate((i_s_edges_coil_S2,i_s_edges_capa_S2,i_s_additionnal_S2),axis=0)
        j_s_S2 = np.concatenate((j_s_edges_coil_S2,j_s_edges_capa_S2,j_s_additionnal_S2),axis=0)
        data_S2 = np.concatenate((data_edges_coil_S2,data_edges_capa_S2,data_additionnal_S2),axis=0)

        mask_res_mutual = ~np.isin(i_s_additionnal_S1,i_s_additionnal_S2)

        i_s_init = np.concatenate((i_s_nodes_S1,i_s_edges_coil_S2,i_s_edges_res_S1,i_s_edges_capa_S2,i_s_mass_S1),axis=0)
        j_s_init = np.concatenate((j_s_nodes_S1,j_s_edges_coil_S2,j_s_edges_res_S1,j_s_edges_capa_S2,j_s_mass_S1),axis=0)
        data_init = np.concatenate((data_nodes_S1,data_edges_coil_S2,data_edges_res_S1,data_edges_capa_S2,data_mass_S1),axis=0)

        self.S_init=(data_init,(i_s_init,j_s_init))

        self.S1 = (data_S1,(i_s_S1.astype(int),j_s_S1.astype(int)))
        self.S2 = (data_S2,(i_s_S2.astype(int),j_s_S2.astype(int)))
        return self.S1,self.S2,self.S_init


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
            node where the tension is enforced
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
        (data_init,(i_s_init,j_s_init)) = self.S_init
        (data_S1,(i_s_S1,j_s_S1)) = self.S1
        (data_S2,(i_s_S2,j_s_S2)) = self.S2
        (data_rhs,(nodes,)) = self.rhs
        new_eq = self.number_intensities+self.size +self.source_count
        sign = np.sign(output_node-input_node)

        ## Adding the new equation in S_init
        if input_node in self.deleted_equation_current:
            data_init = np.append(data_init,np.array([1.,-1.,-sign]),axis=0)
            i_s_init = np.append(i_s_init,np.array([new_eq,new_eq,output_node+self.rescaler[output_node]]),axis=0)
            j_s_init = np.append(j_s_init,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq]),axis=0)
        elif output_node in self.deleted_equation_current:
            data_init = np.append(data_init,np.array([1.,-1.,sign]),axis=0)
            i_s_init = np.append(i_s_init,np.array([new_eq,new_eq,input_node+self.rescaler[input_node]]),axis=0)
            j_s_init = np.append(j_s_init,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq]),axis=0)
        else:
            data_init = np.append(data_init,np.array([1.,-1.,sign,-sign]),axis=0)
            i_s_init = np.append(i_s_init,np.array([new_eq,new_eq,input_node+self.rescaler[input_node],output_node+self.rescaler[output_node]]),axis=0)
            j_s_init = np.append(j_s_init,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq,new_eq]),axis=0)

        ## Adding the new equation in S1
        if input_node in self.deleted_equation_current:
            data_S1 = np.append(data_S1,np.array([1.,-1.,-sign]),axis=0)
            i_s_S1 = np.append(i_s_S1,np.array([new_eq,new_eq,output_node+self.rescaler[output_node]]),axis=0)
            j_s_S1 = np.append(j_s_S1,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq]),axis=0)
        elif output_node in self.deleted_equation_current:
            data_S1 = np.append(data_S1,np.array([1.,-1.,sign]),axis=0)
            i_s_S1 = np.append(i_s_S1,np.array([new_eq,new_eq,input_node+self.rescaler[input_node]]),axis=0)
            j_s_S1 = np.append(j_s_S1,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq]),axis=0)
        else:
            data_S1 = np.append(data_S1,np.array([1.,-1.,sign,-sign]),axis=0)
            i_s_init = np.append(i_s_init,np.array([new_eq,new_eq,input_node+self.rescaler[input_node],output_node+self.rescaler[output_node]]),axis=0)
            j_s_S1 = np.append(j_s_S1,np.array([self.offset_j + input_node,self.offset_j+output_node,new_eq,new_eq]),axis=0)


        ## second member value
        nodes= np.append(nodes, [new_eq],axis=0)
        data_rhs= np.append(data_rhs, [tension],axis=0)


        ## reaffecting the systems and second member
        self.S_init = (data_init,(i_s_init,j_s_init))
        self.S1 = (data_S1,(i_s_S1,j_s_S1))
        self.S2 = (data_S2,(i_s_S2,j_s_S2))
        self.rhs = (data_rhs,(nodes,))

        ## counter for tension source
        self.source_count +=1
        self.source_signs=np.append(self.source_signs,[sign])
        return self.rhs

    def get_init_system(self):
        sys = coo_matrix(self.S_init,shape=(self.number_intensities+self.size+self.source_count,self.number_intensities+self.size+self.source_count))
        rhs = np.zeros(self.number_intensities+self.size+self.source_count)
        (data_rhs,(nodes,)) = self.rhs
        rhs[nodes]=data_rhs
        return sys,rhs

    def get_system(self):
        sys1 = coo_matrix(self.S1,shape=(self.number_intensities+self.size+self.source_count,self.number_intensities+self.size+self.source_count))
        sys2 = coo_matrix(self.S2,shape=(self.number_intensities+self.size+self.source_count,self.number_intensities+self.size+self.source_count))
        rhs = np.zeros(self.number_intensities+self.size+self.source_count)
        (data_rhs,(nodes,)) = self.rhs
        rhs[nodes]=data_rhs
        return sys1,sys2,rhs

    def get_frequency_system(self,omega):
        sys1 = coo_matrix(self.S1,shape=(self.number_intensities+self.size+self.source_count,self.number_intensities+self.size+self.source_count))
        sys2 = coo_matrix(self.S2,shape=(self.number_intensities+self.size+self.source_count,self.number_intensities+self.size+self.source_count))
        rhs = np.zeros(self.number_intensities+self.size+self.source_count)
        (data_rhs,(nodes,)) = self.rhs
        rhs[nodes]=data_rhs
        return sys1+1j*omega*sys2,rhs


    def build_intensity_and_voltage_from_vector(self,sol):
        sign = np.sign(self.all_coords[1]-self.all_coords[0])
        offset_coil = self.coil_data.shape[0]
        offset_res = self.res_data.shape[0]
        offset_capa = self.capa_data.shape[0]
        if self.source_count!=0:
            return (sol[:offset_coil]*sign[:offset_coil],
                sol[offset_coil:offset_coil+offset_res]*sign[offset_coil:offset_coil+offset_res],
                sol[offset_coil+offset_res:offset_coil+offset_res+offset_capa]*sign[offset_coil+offset_res:offset_coil+offset_res+offset_capa],
                sol[self.number_intensities:-self.source_count],
                sol[-self.source_count:]*self.source_signs
            )
        else:
            return (sol[:offset_coil]*sign[:offset_coil],
                sol[offset_coil:offset_coil+offset_res]*sign[offset_coil:offset_coil+offset_res],
                sol[offset_coil+offset_res:offset_coil+offset_res+offset_capa]*sign[offset_coil+offset_res:offset_coil+offset_res+offset_capa],
                sol[self.number_intensities:-self.source_count],
                np.array([],dtype=float)
            )

