from ElecSolver import NetlistParser
import unittest

class TestNetlistParser(unittest.TestCase):
    def setUp(self):
        # Create a temporary netlist file for testing
        self.test_file_path = 'test.net'
        self.wrong_file_type_path = 'test.txt'
        self.param_netlist = "param.cir"
        self.validation_netlist = "validation.cir"
        self.source_netlist = "sources.cir"
        self.precision_parse=12
        with open(self.test_file_path, 'w') as f:
            f.write("R1 n1 n2 10k\n")
            f.write("r2 0 n2 1k\n")
            f.write("L1 n32 n4 100u\n")
            f.write("c1 n5 n6 1n\n")
            f.write("C2 28 q6 1n\n")
            f.write("K1 L1 l2 0.5\n")
            f.write("k12 L122 L21 -0.2\n")
            f.write("l2 n23 n42 10u\n")

        with open(self.wrong_file_type_path, 'w') as f:
            f.write("This is not a netlist file.\n")
            f.write("R1 n1 n2 10k\n")

        with open(self.param_netlist, 'w') as f:
            f.write("R1 n1 n2 {Rt}\n")
            f.write("r2 0 n2 {Rt}\n")
            f.write("L1 n32 n4 100u\n")
            f.write("c1 n5 n6 {Capa}\n")
            f.write("C2 28 q6 {Capa}\n")
            f.write("K1 L1 l2 0.5\n")
            f.write("l2 n23 n42 10u\n")
            f.write(".param Capa=1n\n")
            f.write(".param Rt = 1e3")

        with open(self.validation_netlist, 'w') as f:
            f.write(".param FREQ = 9k\n")
            f.write(".param DUTY = 0.5\n")
            f.write(".param V_LOW = 0\n")
            f.write(".param V_HIGH = 1\n")
            f.write(".param T_PER = {1/FREQ}\n")
            f.write(".param T_ON = {T_PER * DUTY}\n")
            f.write(".param T_RISE = 10n\n")
            f.write(".param T_FALL = 10n\n")
            f.write("\n")
            f.write("V1 in 0 PULSE({V_LOW} {V_HIGH} 0 {T_RISE} {T_FALL} {T_ON} {T_PER})\n")
            f.write("Lfilter in n1 0.5e-3\n")
            f.write("Cfilter n1 0 50e-6\n")
            f.write("Rload n1 0 50\n")
            f.write("Bload n1 0 I=V(Rload)/50\n")

        with open(self.source_netlist, 'w') as f:
            f.write("I1 input 0 DC 1\n")
            f.write("V1 input output 5\n")
            f.write("BI1 output 0 I = V(V1)\n")
            f.write("BV1 sense 0 v=V(output)\n")

    def tearDown(self):
        # Remove the test file after tests
        import os
        if os.path.isfile(self.test_file_path):
            os.remove(self.test_file_path)
        if os.path.isfile(self.wrong_file_type_path):
            os.remove(self.wrong_file_type_path)
        if os.path.isfile(self.param_netlist):
            os.remove(self.param_netlist)
        if os.path.isfile(self.validation_netlist):
            os.remove(self.validation_netlist)
        if os.path.isfile(self.source_netlist):
            os.remove(self.source_netlist)

    def test_init_invalid_file(self):
        """Test initialization with an invalid file path."""
        with self.assertRaises(FileNotFoundError):
            NetlistParser('invalid_path.net')

    def test_init_invalid_file_type(self):
        """Test initialization with an invalid file type."""
        with self.assertRaises(ValueError):
            NetlistParser(self.wrong_file_type_path)

    def test_parse_netlist(self):
        """Test parsing of the netlist file."""
        parser = NetlistParser(self.test_file_path)
        parser._parse_netlist()
        dict_test = {'R1': {'n1': 'n1', 'n2': 'n2', 'value': '10k'},
                      'r2': {'n1': '0', 'n2': 'n2', 'value': '1k'},
                      'L1': {'n1': 'n32', 'n2': 'n4', 'value': '100u'},
                      'c1': {'n1': 'n5', 'n2': 'n6', 'value': '1n'},
                      'C2': {'n1': '28', 'n2': 'q6', 'value': '1n'},
                      'K1': {'L1': 'L1', 'L2': 'l2', 'k': '0.5'},
                      'k12': {'L1': 'L122', 'L2': 'L21', 'k': '-0.2'},
                      'l2': {'n1': 'n23', 'n2': 'n42', 'value': '10u'}}

        for key, value in dict_test.items():
            with self.subTest(key=key):
                if key.startswith(('R', 'r')):
                    self.assertIn(key, parser.resistors)
                    self.assertEqual(parser.resistors[key], value)
                elif key.startswith(('L', 'l')):
                    self.assertIn(key, parser.inductors)
                    self.assertEqual(parser.inductors[key], value)
                elif key.startswith(('C', 'c')):
                    self.assertIn(key, parser.capacitors)
                    self.assertEqual(parser.capacitors[key], value)
                elif key.startswith(('K', 'k')):
                    self.assertIn(key, parser.couplings)
                    self.assertEqual(parser.couplings[key], value)


    def test_parse_si_value(self):
        """Test parsing of SI values."""
        parser = NetlistParser(self.test_file_path)
        dict_test = {'1n': 1e-9, '100u': 100e-6, '10.0µ': 10e-6,
            '2.5m': 2.5e-3, '0.5': 0.5, '10k': 10.e3, '1K': 1.e3,
            '10.5Meg': 10.5e6, '10.5meg': 10.5e6, '10.5M': 10.5e6,
            '4.5G': 4.5e9, '9.5T': 9.5e12,'1.5e3': 1.5e3, '1.5e-3': 1.5e-3,
            '1.5e+3': 1.5e3}
        # Test valid SI values
        for key, value in dict_test.items():
            with self.subTest(key=key):
                self.assertAlmostEqual(parser._parse_si_value(key), value, places=self.precision_parse)


    def test_parse_si_value_invalid(self):
        """Test whether parsing invalid SI values raises errors."""
        # Test invalid SI values
        parser = NetlistParser(self.test_file_path)
        list_wrong_val = ['invalid', '10.5X', '1.5e', '1.5e+', '1.5e-']
        for val in list_wrong_val:
            with self.subTest(val=val):
                with self.assertRaises(ValueError):
                    parser._parse_si_value(val)


    def test_node_map(self):
        """Test mapping of nodes to integers."""
        parser = NetlistParser(self.test_file_path)
        parser._parse_netlist()
        parser._map_nodes()
        dict_test = {'0': 0, 'n1': 1, 'n2': 2, 'n32': 3, 'n4': 4,'n23': 5,
                      'n42': 6, 'n5': 7, 'n6': 8, '28': 9, 'q6': 10}
        for key, value in dict_test.items():
            with self.subTest(key=key):
                self.assertEqual(parser.node_map[key], value)

    def test_coupling_map(self):
        """Test mapping of couplings to integers."""
        parser = NetlistParser(self.test_file_path)
        parser._parse_netlist()
        parser._map_nodes()
        parser._map_couplings()
        dict_test = {'K1': {"L_coords" : [0,1],  "value": 1.5811388300841894e-05}}
        for key, value in dict_test.items():
            with self.subTest(key=key):
                self.assertEqual(parser.coupling_map[key]["L_coords"], value["L_coords"])
                self.assertEqual(parser.coupling_map[key]["value"], value["value"])

    def test_parse_param_values(self):
        parser = NetlistParser(self.param_netlist)
        parser._parse_netlist()
        dict_test = {'Capa': 1e-9,
                     'Rt': 1e3}
        for key, value in dict_test.items():
            with self.subTest(key=key):
                self.assertIn(key, parser.params)
                self.assertEqual(parser.params[key], value)

    def test_composed_param_values_and_sources(self):
        parser = NetlistParser(self.validation_netlist)
        parser.map_netlist()

        expected_params = {
            'FREQ': 9e3,
            'DUTY': 0.5,
            'V_LOW': 0.0,
            'V_HIGH': 1.0,
            'T_PER': 1 / 9e3,
            'T_ON': 0.5 / 9e3,
            'T_RISE': 10e-9,
            'T_FALL': 10e-9,
        }
        for name, value in expected_params.items():
            with self.subTest(name=name):
                self.assertAlmostEqual(
                    parser.params[name], value, places=self.precision_parse)

        self.assertEqual(parser.dipole_map['Lfilter']['value'], 0.5e-3)
        self.assertEqual(parser.dipole_map['Cfilter']['value'], 50e-6)
        self.assertEqual(parser.dipole_map['Rload']['value'], 50.0)
        self.assertEqual(parser.voltage, {'V1': {'n1': 'in', 'n2': '0'}})
        self.assertEqual(parser.current, {'Bload': {'n1': 'n1', 'n2': '0'}})
        self.assertEqual(parser.dipole_map['V1'], {'nodes': [2, 0]})
        self.assertEqual(parser.dipole_map['Bload'], {'nodes': [1, 0]})

    def test_independent_and_behavioral_sources_share_attributes(self):
        parser = NetlistParser(self.source_netlist)
        parser.map_netlist()

        self.assertEqual(
            parser.current,
            {
                'I1': {'n1': 'input', 'n2': '0'},
                'BI1': {'n1': 'output', 'n2': '0'},
            },
        )
        self.assertEqual(
            parser.voltage,
            {
                'V1': {'n1': 'input', 'n2': 'output'},
                'BV1': {'n1': 'sense', 'n2': '0'},
            },
        )
        for source in ('I1', 'BI1', 'V1', 'BV1'):
            with self.subTest(source=source):
                self.assertIn(source, parser.dipole_map)
                self.assertNotIn('value', parser.dipole_map[source])

    def test_param_expression_rejects_forward_reference(self):
        with open(self.test_file_path, 'w') as f:
            f.write(".param FIRST = {SECOND * 2}\n")
            f.write(".param SECOND = 10\n")

        parser = NetlistParser(self.test_file_path)
        with self.assertRaisesRegex(
                ValueError, "Unknown or forward parameter reference 'SECOND'"):
            parser._parse_netlist()

    def test_composed_param_values_map_to_passive_components(self):
        with open(self.test_file_path, 'w') as f:
            f.write(".param BASE = 10\n")
            f.write(".param SCALE = 2\n")
            f.write("Rexpr n1 0 {BASE * SCALE}\n")
            f.write("Lexpr n1 0 {BASE / SCALE}\n")
            f.write("Cexpr n1 0 {(BASE + SCALE) * 1n}\n")

        parser = NetlistParser(self.test_file_path)
        parser.map_netlist()

        self.assertEqual(parser.dipole_map['Rexpr']['value'], 20.0)
        self.assertEqual(parser.dipole_map['Lexpr']['value'], 5.0)
        self.assertAlmostEqual(
            parser.dipole_map['Cexpr']['value'],
            12e-9,
            places=self.precision_parse,
        )

    def test_param_affectation(self):
        parser = NetlistParser(self.param_netlist)
        parser.map_netlist()
        for dipole in parser.dipole_map:
            value = parser.dipole_map[dipole]["value"]
            if dipole.startswith(("R","r")):
                self.assertEqual(value, 1000.)
            if dipole.startswith(("C","c")):
                self.assertEqual(value, 1.e-9)
    
    def test_generate_system(self):
        parser = NetlistParser(self.param_netlist)
        parser.map_netlist()
        parser.generate_temporal_system()



if __name__ == '__main__':
    unittest.main()
