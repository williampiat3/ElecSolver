import numpy as np
from .NetlistGenerator import ResistorNetlist, CapacitorNetlist, InductorNetlist, MutualInductance, SubCircuitDefinition
from .TemporalSystemBuilder import TemporalSystemBuilder

class TemporalNetlistDumper(object): 
    def __init__(self, temporal_system_builder: TemporalSystemBuilder):
        """Constructor of TemporalNetlistDumper. Allows dumping TemporalSystemBuilder as subckt compatible with spice-like tools

        Args:
            temporal_system_builder (TemporalSystemBuilder): TemporalSystemBuilder from elecsolv to be exported
        """
        if not isinstance(temporal_system_builder, TemporalSystemBuilder):
            raise TypeError("temporal_system_builder must be an instance of "+ TemporalSystemBuilder.__name__)
        self.system_builder = temporal_system_builder
        self.ports_map ={} 
    
    def add_port(self, port_number: int, port_name: str):
        """Add port to subckt definition, will be exposed and accessible in spice tools

        Args:
            port_number (int): port number used in temporal_system_builder to be used as port
            port_name (str): port name to be used for export 
        """ 
        self.ports_map.update({port_number:port_name})

    def generate_subcircuit_file(self, subcircuit_name : str = "esolvsub", file_name : str ="esolvsub.net"):
        """Dump a text file with the the definition of the subckt with provided name and file name

        Args:
            subcircuit_name (str, optional): name of the subcircuit to be generated. Defaults to "esolvsub".
            file_name (str, optional): file name to dump the generate the netlist to. Defaults to "esolvsub.net".
        """
        if not bool(self.ports_map):  # check if dict is empty 
            raise ValueError("ports_map is empty, cannot generate subckt")
        list_res =[]
        list_coil=[]
        list_capa=[]
        list_k=[]
        list_out = [list_res,list_coil,list_capa]
        list_net = [ResistorNetlist,InductorNetlist,CapacitorNetlist]
        list_data = [self.system_builder.res_data,self.system_builder.coil_data, self.system_builder.capa_data]
        list_coords = [self.system_builder.res_coords,self.system_builder.coil_coords, self.system_builder.capa_coords]

        for index_type_dipole, current_class in enumerate(list_net):
            curr_list_out = list_out[index_type_dipole]
            curr_data = list_data[index_type_dipole]
            curr_coords = list_coords[index_type_dipole]

            for index, value in  enumerate(curr_data):
                node1 = curr_coords[0][index]
                node2 = curr_coords[1][index]

                if node1 in self.ports_map:
                    node1=self.ports_map[node1]
                else:
                    node1=str(node1+1) # +1 because 0 node is gnd in spice !

                if node2 in self.ports_map:
                    node2=self.ports_map[node2]
                else:
                    node2=str(node2+1) # +1 because 0 node is gnd in spice !
                curr_list_out.append(current_class(str(index), node1, node2,  value))
        
        for index_mutual, mutual in enumerate(self.system_builder.inductive_mutuals_data):
            l1 = self.system_builder.inductive_mutuals_coords[0][index_mutual]
            l2 = self.system_builder.inductive_mutuals_coords[1][index_mutual]
            l1_val = self.system_builder.coil_data[l1]
            l2_val = self.system_builder.coil_data[l2]
            k = mutual/np.sqrt(l1_val*l2_val)
            if k!=0:
                list_k.append(MutualInductance(str(index_mutual),"L"+str(l1), "L"+str(l2),k))

        flat_list=[]
        for row in list_out:
            flat_list += row
        flat_list+=list_k
        subskt = SubCircuitDefinition(subcircuit_name, flat_list, self.ports_map.values())
        with open(file_name, "w", encoding="utf-8") as f:
            f.write("*Generated subckt from ElecSolver\n")
            f.write(subskt._char_netlist())


