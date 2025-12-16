import re
import os
from numpy import sqrt, array, concatenate, arange
from .TemporalSystemBuilder import TemporalSystemBuilder

class NetlistParser():
    """
    A class to parse a netlist file and extract component values and connections.

    Attributes
    ----------
    file_path : str
        Path to the netlist file.
    node_map : dict
        Map of node names to unique integers.
    resistors : dict
        Dictionary of resistors with their connections and values.
    inductors : dict
        Dictionary of inductors with their connections and values.
    capacitors : dict
        Dictionary of capacitors with their connections and values.
    couplings : dict
        Dictionary of couplings with their connections and values.
    """
    # Regex patterns
    RESISTOR_PATTERN  = r'^([Rr][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    INDUCTOR_PATTERN  = r'^([Ll][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    CAPACITOR_PATTERN = r'^([Cc][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    COUPLING_PATTERN  = r'^([Kk][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    REALCOUPLING_PATTERN  = r'^([Ww][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    PARAM_PATTERN = r'^\.param\s+([\w]+)\s*=\s*(.+)$'

    SI_COEF = {
        'f': 1e-15,
        'p': 1e-12,
        'n': 1e-9,
        'u': 1e-6,
        'µ': 1e-6,  # Allow unicode micro
        'm': 1e-3,
        '': 1,
        'k': 1e3,
        'K': 1e3,
        'meg': 1e6,
        'Meg': 1e6,
        'M': 1e6,
        'g': 1e9,
        'G': 1e9,
        't': 1e12,
        'T': 1e12
    }

    def __init__(self, file_path):
        """
        Initialize the NetlistParser with a file path.
        Checks whether the file exists and is a valid netlist file.

        Parameters
        ----------
        file_path : str
            Path to the netlist file.
        """
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        if not file_path.endswith(('.net', '.cir', '.sp')):
            raise ValueError(f"Invalid file type: {file_path}. Expected a .net file.")

        self.file_path = file_path
        self.node_map = {}
        self.dipole_map = {}
        self.coupling_map = {}
        self.real_coupling_map = {}
        self.params = {}
        self.max_index_node = 0


    def map_netlist(self):
        """
        Parse the netlist file and create mappings for components and nodes.
        This method will parse the netlist file, extract component values and connections,
        and create a mapping of node names to unique integers.
        """
        self._parse_netlist()
        self._map_nodes()
        self._map_couplings()
        #print(self.dipole_map, self.coupling_map)


    def _parse_netlist(self):
        """
        Parse the netlist file and extract component values and connections.
        """

        with open(self.file_path, 'r') as file:
            data = file.read()

        # Use re.MULTILINE to process line by line without looping
        self.resistors = {m[0]: {'n1': m[1], 'n2': m[2], 'value': m[3]}
                 for m in re.findall(self.RESISTOR_PATTERN, data, re.MULTILINE)}
        self.inductors = {m[0]: {'n1': m[1], 'n2': m[2], 'value': m[3]}
                 for m in re.findall(self.INDUCTOR_PATTERN, data, re.MULTILINE)}
        self.capacitors = {m[0]: {'n1': m[1], 'n2': m[2], 'value': m[3]}
                  for m in re.findall(self.CAPACITOR_PATTERN, data, re.MULTILINE)}
        self.couplings = {m[0]: {'L1': m[1], 'L2': m[2], 'k': m[3]}
                 for m in re.findall(self.COUPLING_PATTERN, data, re.MULTILINE)}
        self.real_couplings = {m[0]: {'L1': m[1], 'L2': m[2], 'k': m[3]}
                 for m in re.findall(self.REALCOUPLING_PATTERN, data, re.MULTILINE)}
        self.params = {m[0]: m[1]
                for m in re.findall(self.PARAM_PATTERN, data, re.MULTILINE)}



    def _parse_si_value(self, value_str):
        """
        Parse a string representing a value with an SI prefix and convert it to a float.

        Parameters
        ----------
        value_str : str
            The string to parse.

        Returns
        -------
        float
            The parsed value in standard (base-unit) form.

        Raises
        ------
        ValueError
            If the string is not a valid SI value.
        """
        try:
            return float(value_str)
        except ValueError:
            pass
        if value_str.startswith("{") and value_str.endswith("}"):
            return self._parse_si_value(self.params[value_str[1:-1]])

        match = re.fullmatch(r'\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)([a-zA-Zµ]*)\s*', value_str)
        if not match:
            raise ValueError(f"Invalid SI value format: '{value_str}'")

        number, prefix = match.groups()

        # First try exact match (e.g. 'M', 'Meg')
        if prefix in self.SI_COEF:
            multiplier = self.SI_COEF[prefix]
        # Then try lowercase for variants like 'meg'
        elif prefix.lower() in self.SI_COEF:
            multiplier = self.SI_COEF[prefix.lower()]
        else:
            raise ValueError(f"Unknown SI prefix: '{prefix}'")

        return float(number) * multiplier

    def _map_nodes(self):
        """
        Create a mapping of node names to unique integers.
        This is used for circuit simulation and analysis.
        """
        # Create a mapping of node names to unique integers
        self.node_map = {'0': 0}  # Ground node
        self.dipole_map = {}
        dipole_list = [self.resistors, self.inductors, self.capacitors]
        index_node=1
        for component in dipole_list:
            for dipole in component:
                n1 = component[dipole]['n1']
                n2 = component[dipole]['n2']
                if n1 not in self.node_map:
                    self.node_map[n1] = index_node
                    index_node += 1
                if n2 not in self.node_map:
                    self.node_map[n2] = index_node
                    index_node += 1
        self.max_index_node = index_node
        # Create a mapping of dipole names to node pairs
        for component in dipole_list:
            for dipole in component:
                self.dipole_map[dipole] = {"nodes":[self.node_map[component[dipole]['n1']],
                                          self.node_map[component[dipole]['n2']]],
                                          "value":self._parse_si_value(component[dipole]['value'])}

    def _map_couplings(self):
        """
        Create a mapping of coupling names to unique integers nodes.
        This is used for circuit simulation and analysis.
        """
        self.coupling_map = {}
        list_l_name = list(self.inductors.keys())
        for coupling in self.couplings:
            L1 = self.couplings[coupling]['L1']
            L2 = self.couplings[coupling]['L2']
            if L1 not in self.inductors:
                pass
            elif L2 not in self.inductors:
                pass
            else:
                M=self._parse_si_value(self.couplings[coupling]['k'])*sqrt(self.dipole_map[L1]['value']*self.dipole_map[L2]['value'])
                self.coupling_map[coupling]={"L_coords":[list_l_name.index(L1), list_l_name.index(L2)],"value":M}

        self.real_coupling_map = {}
        for coupling in self.real_couplings:
            L1 = self.real_couplings[coupling]['L1']
            L2 = self.real_couplings[coupling]['L2']
            if L1 not in self.inductors:
                pass
            elif L2 not in self.inductors:
                pass
            else:
                Rij=self._parse_si_value(self.real_couplings[coupling]['k'])
                self.real_coupling_map[coupling]={"L_coords":[list_l_name.index(L1), list_l_name.index(L2)],"value":Rij}

    def _fill_array_circuit(self,indexes, values, dipole_list):
        """
        Convert dipole mapping into a format compatible with elecsolve.

        Parameters
        ----------
        indexes : list of int
            Coordinates for elecsolve nodes to be filled (length 2 * n).
        values : list of float
            Values corresponding to the coordinates.
        dipole_list : list of str
            Names of dipoles to be added.

        Returns
        -------
        indexes : list of int
            Filled coordinates for elecsolve nodes (length 2 * n).
        values : list of float
            Filled values matching the returned coordinates.
        """
        for dipole in dipole_list:
            nodes = self.dipole_map[dipole]["nodes"]
            value = self.dipole_map[dipole]["value"]
            indexes = concatenate((indexes, [[nodes[0]],[nodes[1]]]), axis=1)
            values = concatenate((values, [value]), axis=0)
        return indexes, values

    def _fill_array_coupling(self,indexes, values, coupling_map):
        """
        Convert dipole coupling into a format compatible with elecsolve.

        Parameters
        ----------
        indexes : list of int
            Coordinates for elecsolve nodes to be filled (length 2 * n).
        values : list of float
            Values corresponding to the coordinates.
        coupling_map : list of str
            Names of couplings to be added.

        Returns
        -------
        indexes : list of int
            Filled coordinates for elecsolve nodes (length 2 * n).
        values : list of float
            Filled values matching the returned coordinates.
        """
        for coupling in coupling_map:
            nodes = coupling_map[coupling]["L_coords"]
            value = coupling_map[coupling]["value"]
            indexes = concatenate((indexes, [[nodes[0]],[nodes[1]]]), axis=1)
            values = concatenate((values, [value]), axis=0)
        return indexes, values

    def generate_temporal_system(self):
        """
        Transform a netlist into a matrix representation for use with the Python sparse solver.

        Returns
        -------
        elec_sys : TemporalSystemBuilder
            TemporalSystemBuilder instance corresponding to the provided netlist.
        """
        index_res, resistors = array([[],[]], dtype=int), array([],dtype=float)
        index_coils, coils_data = array([[],[]], dtype=int), array([],dtype=float)
        capa_coords, capa_data = array([[],[]], dtype=int), array([],dtype=float)
        coords_mutual, data_mutual = array([[],[]], dtype=int), array([],dtype=float)
        res_mutuals_coords,res_mutuals_data = array([[],[]], dtype=int), array([],dtype=float)

        index_res, resistors = self._fill_array_circuit(index_res, resistors, self.resistors)
        index_coils, coils_data = self._fill_array_circuit( index_coils, coils_data, self.inductors)
        capa_coords, capa_data = self._fill_array_circuit(capa_coords, capa_data, self.capacitors)

        coords_mutual, data_mutual = self._fill_array_coupling(coords_mutual, data_mutual , self.coupling_map)
        res_mutuals_coords,res_mutuals_data = self._fill_array_coupling(res_mutuals_coords,res_mutuals_data, self.real_coupling_map)

        elec_sys = TemporalSystemBuilder(index_coils,coils_data,index_res,resistors,capa_coords,capa_data,coords_mutual,data_mutual,res_mutuals_coords,res_mutuals_data)

        return elec_sys