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
