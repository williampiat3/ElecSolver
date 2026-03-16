from ElecSolver import TemporalSystemBuilder 
from ElecSolver import TemporalNetlistDumper
import numpy as np
import unittest

class TestNetlistDumper(unittest.TestCase):
    def setUp(self):
        # Create a temporary netlist file for testing
        self.test_file_path = 'test_gen.net'

        # build simple circuit
        ## Dipole build
        res_coords  = np.array([[0,2],[1,3]],dtype=int)
        res_data = np.array([1,1],dtype=float)
        coil_coords  = np.array([[1,0],[2,3]],dtype=int)
        coil_data = np.array([1,1],dtype=float)
        capa_coords = np.array([[1,2],[3,0]],dtype=int)
        capa_data = np.array([1,1],dtype=float)
        ## coupling factors
        mutuals_coords=np.array([[0],[1]],dtype=int)
        mutuals_data = np.array([1e-5],dtype=float)
        res_mutuals_coords=np.array([[],[]],dtype=int)
        res_mutuals_data = np.array([],dtype=float)

        self.elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)

    def tearDown(self):
        # Remove the test file after tests
        import os
        if os.path.isfile(self.test_file_path):
            os.remove(self.test_file_path)

    def test_init_invalid_system(self):
        with self.assertRaises(TypeError):
            TemporalNetlistDumper("self.elec_sys")

    def test_init_valid_system(self):
        TemporalNetlistDumper(self.elec_sys)
    
    def test_add_port(self):
        test = TemporalNetlistDumper(self.elec_sys)
        test.add_port(1,'in')
        self.assertDictEqual(test.ports_map, {1:'in'})

    def test_fail_dump_no_port_defined(self):
        test = TemporalNetlistDumper(self.elec_sys)
        with self.assertRaises(ValueError):
            test.generate_subcircuit_file()

    def  test_dump_with_port_defined(self):
        test = TemporalNetlistDumper(self.elec_sys)
        test.add_port(1,'in')
        test.add_port(0,'out')
        test.generate_subcircuit_file(file_name= self.test_file_path)

    
    


    

    





#   def test_init_invalid_file_type(self):
#       """Test initialization with an invalid file type."""
#       with self.assertRaises(ValueError):
#           NetlistParser(self.wrong_file_type_path)
#
#   def test_parse_netlist(self):
#       """Test parsing of the netlist file."""
#       parser = NetlistParser(self.test_file_path)
#       parser._parse_netlist()
#       dict_test = {'R1': {'n1': 'n1', 'n2': 'n2', 'value': '10k'},
#                     'r2': {'n1': '0', 'n2': 'n2', 'value': '1k'},
#                     'L1': {'n1': 'n32', 'n2': 'n4', 'value': '100u'},
#                     'c1': {'n1': 'n5', 'n2': 'n6', 'value': '1n'},
#                     'C2': {'n1': '28', 'n2': 'q6', 'value': '1n'},
#                     'K1': {'L1': 'L1', 'L2': 'l2', 'k': '0.5'},
#                     'k12': {'L1': 'L122', 'L2': 'L21', 'k': '-0.2'},
#                     'l2': {'n1': 'n23', 'n2': 'n42', 'value': '10u'}}
#
#       for key, value in dict_test.items():
#           with self.subTest(key=key):
#               if key.startswith(('R', 'r')):
#                   self.assertIn(key, parser.resistors)
#                   self.assertEqual(parser.resistors[key], value)
#               elif key.startswith(('L', 'l')):
#                   self.assertIn(key, parser.inductors)
#                   self.assertEqual(parser.inductors[key], value)
#               elif key.startswith(('C', 'c')):
#                   self.assertIn(key, parser.capacitors)
#                   self.assertEqual(parser.capacitors[key], value)
#               elif key.startswith(('K', 'k')):
#                   self.assertIn(key, parser.couplings)
#                   self.assertEqual(parser.couplings[key], value)
#
#
#   def test_parse_si_value(self):
#       """Test parsing of SI values."""
#       parser = NetlistParser(self.test_file_path)
#       dict_test = {'1n': 1e-9, '100u': 100e-6, '10.0µ': 10e-6,
#           '2.5m': 2.5e-3, '0.5': 0.5, '10k': 10.e3, '1K': 1.e3,
#           '10.5Meg': 10.5e6, '10.5meg': 10.5e6, '10.5M': 10.5e6,
#           '4.5G': 4.5e9, '9.5T': 9.5e12,'1.5e3': 1.5e3, '1.5e-3': 1.5e-3,
#           '1.5e+3': 1.5e3}
#       # Test valid SI values
#       for key, value in dict_test.items():
#           with self.subTest(key=key):
#               self.assertAlmostEqual(parser._parse_si_value(key), value, places=self.precision_parse)
#
#
#   def test_parse_si_value_invalid(self):
#       """Test if parsing of invalid SI values raise errors."""
#       # Test invalid SI values
#       parser = NetlistParser(self.test_file_path)
#       list_wrong_val = ['invalid', '10.5X', '1.5e', '1.5e+', '1.5e-']
#       for val in list_wrong_val:
#           with self.subTest(val=val):
#               with self.assertRaises(ValueError):
#                   parser._parse_si_value(val)
#
#
#   def test_node_map(self):
#       """Test mapping of nodes to integers."""
#       parser = NetlistParser(self.test_file_path)
#       parser._parse_netlist()
#       parser._map_nodes()
#       dict_test = {'0': 0, 'n1': 1, 'n2': 2, 'n32': 3, 'n4': 4,'n23': 5,
#                     'n42': 6, 'n5': 7, 'n6': 8, '28': 9, 'q6': 10}
#       for key, value in dict_test.items():
#           with self.subTest(key=key):
#               self.assertEqual(parser.node_map[key], value)
#
#   def test_coupling_map(self):
#       """Test mapping of couplings to integers."""
#       parser = NetlistParser(self.test_file_path)
#       parser._parse_netlist()
#       parser._map_nodes()
#       parser._map_couplings()
#       dict_test = {'K1': {"L_coords" : [0,1],  "value": 1.5811388300841894e-05}}
#       for key, value in dict_test.items():
#           with self.subTest(key=key):
#               self.assertEqual(parser.coupling_map[key]["L_coords"], value["L_coords"])
#               self.assertEqual(parser.coupling_map[key]["value"], value["value"])
#
#   def test_parse_param_values(self):
#       parser = NetlistParser(self.param_netlist)
#       parser._parse_netlist()
#       dict_test = { 'Capa':'1n',
#                     'Rt': '1e3'}
#       for key, value in dict_test.items():
#           with self.subTest(key=key):
#               self.assertIn(key, parser.params)
#               self.assertEqual(parser.params[key], value)
#   
#   def test_param_affectation(self):
#       parser = NetlistParser(self.param_netlist)
#       parser.map_netlist()
#       for dipole in parser.dipole_map:
#           value = parser.dipole_map[dipole]["value"]
#           if dipole.startswith(("R","r")):
#               self.assertEqual(value, 1000.)
#           if dipole.startswith(("C","c")):
#               self.assertEqual(value, 1.e-9)
#   
#   def test_generate_system(self):
#       parser = NetlistParser(self.param_netlist)
#       parser.map_netlist()
#       parser.generate_temporal_system()



if __name__ == '__main__':
    unittest.main()


