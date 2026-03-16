import datetime


class Netlist(object):
    """
    Represents a netlist for electric circuit simulation in spice format.

    This class serves as a base class for creating specific netlist instances.
    It includes methods for retrieving the type character, generating a string
    representation of the netlist, and creating a full name for the netlist by
    combining its type character and name.

    Parameters:
    - name (str): Name of the netlist.
    - value: Value associated with the netlist.

    Methods:
    - get_type_char(): Abstract method to be implemented by subclasses. Returns
      a character representing the type of the netlist.
    - _char_netlist(): Abstract method to be implemented by subclasses. Returns
      a string representation of the netlist in the spice format.
    - full_name(): Combines the type character and name to create a full name
      for the netlist.

    Note: Subclasses must implement the abstract methods get_type_char() and
     _char_netlist() to define the specific behavior of the netlist.
    """

    def __init__(self, name, value):
        """
        Constructor a Netlist instance with a name and a value.
        :param value:  representing the value of the netlist.
        :param name (str): representing the name of the netlist
        """ 
        self.value = value
        self.name = name   
    
    def get_type_char(self):
        """Abstract method to be implemented by subclasses.
           Returns a character representing the type of the netlist.
        """
        raise NotImplementedError()
    
    def _char_netlist(self):
        """Abstract method to be implemented by subclasses.
           Returns a string representation of the netlist in the spice format.
        """
        raise NotImplementedError()

    def full_name(self):
        """
        Combines the type character and name to create a full name for the netlist.
        """
        return self.get_type_char() + str(self.name)


class DipoleNetlist(Netlist):
    """
     Represents a dipole netlist for electric circuit simulation in spice format.
     
     This class extends the base class Netlist and adds specific attributes for
     dipole netlists. It includes a method to generate a string representation
     of the dipole netlist in the spice format.
     
     Parameters:
     - name (str): Name of the dipole netlist.
     - n1: Node 1 of the dipole netlist.
     - n2: Node 2 of the dipole netlist.
     - value: Value associated with the dipole netlist.
     
     Methods:
     - _char_netlist(): Overrides the method in the base class. Returns a string
       representation of the dipole netlist in the spice format.
     
     Note: The dipole netlist inherits from the base class Netlist.
    """

    def __init__(self, name, n1, n2, value):
        """
        Constructor for a DipoleNetlist instance with a name, nodes, and a value.
        :param name (str): Name of the dipole netlist.
        :param n1: Node 1 of the dipole netlist.
        :param n2: Node 2 of the dipole netlist.
        :param value: Value associated with the dipole netlist.
        """
        super(DipoleNetlist, self).__init__(name, value)
        self.n1 = n1  # Node 1 of the dipole netlist.
        self.n2 = n2  # Node 2 of the dipole netlist.

    def _char_netlist(self):
        """
        Overrides the method in the base class.
        Returns a string representation of the dipole netlist in the spice format.
        """
        net_char = self.full_name() + " " + str(self.n1) + " " + str(self.n2) + " " + str(self.value)
        return net_char

    
class ResistorNetlist(DipoleNetlist):
    """
     Represents a resistor netlist for electric circuit simulation in spice format.
     
     This class extends the DipoleNetlist class and specializes in representing
     resistors. It overrides the get_type_char method to return the resistor type
     character ('R').
     
     Parameters:
     - name (str): Name of the resistor netlist.
     - n1: Node 1 of the resistor netlist.
     - n2: Node 2 of the resistor netlist.
     - value: Value associated with the resistor netlist.
     
     Methods:
     - get_type_char(): Overrides the method in the base class. Returns the
       resistor type character ('R').
     
     Note: The resistor netlist inherits from the DipoleNetlist class.
    """

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: the resistor type character ('R').
        """
        return 'R'


class CapacitorNetlist(DipoleNetlist):
    """
    Represents a capacitor netlist for electric circuit simulation in spice format.

    This class extends the DipoleNetlist class and specializes in representing
    capacitors. It overrides the get_type_char method to return the capacitor type
    character ('C').

    Parameters:
    - name (str): Name of the capacitor netlist.
    - n1: Node 1 of the capacitor netlist.
    - n2: Node 2 of the capacitor netlist.
    - value: Value associated with the capacitor netlist.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      capacitor type character ('C').

    Note: The capacitor netlist inherits from the DipoleNetlist class.
    """

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The capacitor type character ('C').
        """
        return 'C'


class InductorNetlist(DipoleNetlist):
    """
    Represents an inductor netlist for electric circuit simulation in spice format.

    This class extends the DipoleNetlist class and specializes in representing
    inductors. It overrides the get_type_char method to return the inductor type
    character ('L').

    Parameters:
    - name (str): Name of the inductor netlist.
    - n1: Node 1 of the inductor netlist.
    - n2: Node 2 of the inductor netlist.
    - value: Value associated with the inductor netlist.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      inductor type character ('L').

    Note: The inductor netlist inherits from the DipoleNetlist class.
    """

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The inductor type character ('L').
        """
        return 'L'


class VoltageSourceNetlist(DipoleNetlist):
    """
    Represents a voltage source netlist for electric circuit simulation in spice format.

    This class extends the DipoleNetlist class and specializes in representing
    voltage sources. It overrides the get_type_char method to return the voltage
    source type character ('V').

    Parameters:
    - name (str): Name of the voltage source netlist.
    - n1: Node 1 of the voltage source netlist.
    - n2: Node 2 of the voltage source netlist.
    - value: Value associated with the voltage source netlist.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      voltage source type character ('V').

    Note: The voltage source netlist inherits from the DipoleNetlist class.
    """

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The voltage source type character ('V').
        """
        return 'V'


class CurrentSourceNetlist(DipoleNetlist):
    """
    Represents a current source netlist for electric circuit simulation in spice format.

    This class extends the DipoleNetlist class and specializes in representing
    current sources. It overrides the get_type_char method to return the current
    source type character ('I').

    Parameters:
    - name (str): Name of the current source netlist.
    - n1: Node 1 of the current source netlist.
    - n2: Node 2 of the current source netlist.
    - value: Value associated with the current source netlist.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      current source type character ('I').

    Note: The current source netlist inherits from the DipoleNetlist class.
    """

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The current source type character ('I').
        """
        return 'I'


class VcontrolNetlist(DipoleNetlist):
    """
    Represents a voltage-controlled source netlist for electric circuit simulation in spice format.

    This class extends the DipoleNetlist class and specializes in representing
    voltage-controlled sources. It adds specific attributes for the control
    nodes (nplus, nminus) and overrides the _char_netlist method to generate a
    string representation of the netlist in the spice format.

    Parameters:
    - name (str): Name of the voltage-controlled voltage source netlist.
    - n1: Node 1 of the voltage-controlled voltage source netlist.
    - n2: Node 2 of the voltage-controlled voltage source netlist.
    - nplus: Positive control node of the voltage-controlled voltage source.
    - nminus: Negative control node of the voltage-controlled voltage source.
    - gain: Gain associated with the voltage-controlled voltage source.

    Methods:
    - _char_netlist(): Overrides the method in the base class. Returns a string
      representation of the voltage-controlled voltage source netlist in the spice format.

    Note: The VcontrolNetlist inherits from the DipoleNetlist class.
    """

    def __init__(self, name, n1, n2, nplus, nminus, gain):
        """
        Constructor for a VcontrolNetlist instance with a name, nodes, control
        nodes, and a gain.
        :param name (str): Name of the voltage-controlled voltage source netlist.
        :param n1: Node 1 of the voltage-controlled voltage source netlist.
        :param n2: Node 2 of the voltage-controlled voltage source netlist.
        :param nplus: Positive control node of the voltage-controlled voltage source.
        :param nminus: Negative control node of the voltage-controlled voltage source.
        :param gain: Gain associated with the voltage-controlled voltage source.
        """
        super(VcontrolNetlist, self).__init__(name, n1, n2, gain)
        self.nplus = nplus  # Positive control node of the voltage-controlled voltage source.
        self.nminus = nminus  # Negative control node of the voltage-controlled voltage source.

    def _char_netlist(self):
        """
        Overrides the method in the base class.
        Returns a string representation of the voltage-controlled voltage source
        netlist in the spice format.
        """
        char_element = self.get_type_char()
        net_char = char_element + str(self.name) + " " + str(self.n1) + " " + str(self.n2) + " " + str(self.nplus) + " "
        net_char = net_char + str(self.nminus) + " " + str(self.value)
        return net_char


class VoltVcontrolledNet(VcontrolNetlist):
    """
    Represents a voltage-controlled voltage source netlist for electric circuit simulation in spice format.

    This class extends the VcontrolNetlist class and specializes in representing
    voltage-controlled voltage sources. It overrides the get_type_char method to return
    the voltage-controlled voltage source type character ('E').

    Parameters:
    - name (str): Name of the voltage-controlled voltage source netlist.
    - n1: Node 1 of the voltage-controlled voltage source netlist.
    - n2: Node 2 of the voltage-controlled voltage source netlist.
    - nplus: Positive control node of the voltage-controlled voltage source.
    - nminus: Negative control node of the voltage-controlled voltage source.
    - gain: Gain associated with the voltage-controlled voltage source.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      voltage-controlled voltage source type character ('E').

    Note: The VoltVconroledNet inherits from the VcontrolNetlist class.
    """

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The voltage-controlled voltage source type character ('E').
        """
        return 'E'


class AmpsVcontrolledNet(VcontrolNetlist):
    """
    Represents a current-controlled source netlist for electric circuit simulation in spice format.

    This class extends the VcontrolNetlist class and specializes in representing
    voltage-controlled current sources. It overrides the get_type_char method to return
    the voltage-controlled current source type character ('G').

    Parameters:
    - name (str): Name of the voltage-controlled current source netlist.
    - n1: Node 1 of the voltage-controlled current source netlist.
    - n2: Node 2 of the voltage-controlled current source netlist.
    - nplus: Positive control node of the voltage-controlled current source.
    - nminus: Negative control node of the voltage-controlled current source.
    - gain: Gain associated with the voltage-controlled current source.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      voltage-controlled current source type character ('G').

    Note: The AmpsVconroledNet inherits from the VcontrolNetlist class.
    """

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The voltage-controlled current source type character ('G').
        """
        return 'G'


class IcontrolNetlist(DipoleNetlist):
    """
    Represents a current-controlled source netlist for electric circuit simulation in spice format.

    This class extends the DipoleNetlist class and specializes in representing
    current-controlled sources. It adds specific attributes for the
    controlling voltage source (Vcontrol) and overrides the _char_netlist method
    to generate a string representation of the netlist in the spice format.

    Parameters:
    - name (str): Name of the current-controlled voltage source netlist.
    - n1: Node 1 of the current-controlled voltage source netlist.
    - n2: Node 2 of the current-controlled voltage source netlist.
    - Vcontrol: I(Vcontrol) will be used to control the source.
    - gain: Gain associated with the current-controlled source.

    Methods:
    - _char_netlist(): Overrides the method in the base class. Returns a string
      representation of the current-controlled voltage source netlist in the spice format.

    Note: The IcontrolNetlist inherits from the DipoleNetlist class.
    """

    def __init__(self, name, n1, n2, Vcontrol, gain):
        """
        Constructor for an IcontrolNetlist instance with a name, nodes, a
        controlling voltage source, and a gain.
        :param name (str): Name of the current-controlled voltage source netlist.
        :param n1: Node 1 of the current-controlled voltage source netlist.
        :param n2: Node 2 of the current-controlled voltage source netlist.
        :param Vcontrol: I(Vcontrol) control the value of the current-controlled source.
        :param gain: Gain associated with the current-controlled source.
        """
        super(IcontrolNetlist, self).__init__(name, n1, n2, gain)
        self.Vcontrol = Vcontrol  # Controlling voltage source of the current-controlled voltage source.

    def _char_netlist(self):
        """
        Overrides the method in the base class.
        Returns a string representation of the current-controlled voltage source
        netlist in the spice format.
        """
        char_element = self.get_type_char()
        net_char = char_element + str(self.name) + " " + str(self.n1) + " " + str(self.n2) + " " + str(self.Vcontrol) + " "
        net_char = net_char + " " + str(self.value)
        return net_char


class VoltageCurrentControlled(IcontrolNetlist):
    """
    Represents a current-controlled voltage source netlist for electric circuit simulation in spice format.

    This class extends the IcontrolNetlist class and specializes in representing
    voltage-controlled current sources. It overrides the get_type_char method to return
    the current-controlled voltage source type character ('H').

    Parameters:
    - name (str): Name of the voltage-controlled current source netlist.
    - n1: Node 1 of the voltage-controlled current source netlist.
    - n2: Node 2 of the voltage-controlled current source netlist.
    - Vcontrol: I(Vcontrol) will be used to control the source.
    - gain: Gain associated with the voltage-controlled current source.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      voltage-controlled current source type character ('H').

    Note: The VoltageCurrentControlled inherits from the IcontrolNetlist class.
    """

    def __init__(self, name, n1, n2, Vcontrol, gain):
        """
        Constructor for a VoltageCurrentControlled instance with a name, nodes, a
        controlling voltage source, and a gain.
        :param name (str): Name of the voltage-controlled current source netlist.
        :param n1: Node 1 of the voltage-controlled current source netlist.
        :param n2: Node 2 of the voltage-controlled current source netlist.
        :param Vcontrol: I(Vcontrol) control the value of the current-controlled voltage source
        :param gain: Gain associated with the voltage-controlled current source.
        """
        super(VoltageCurrentControlled, self).__init__(name, n1, n2, Vcontrol, gain)

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The voltage-controlled current source type character ('H').
        """
        return 'H'

    
class CurrentCurrentControlled(IcontrolNetlist):
    """
    Represents a current-controlled current source netlist for electric circuit simulation in spice format.

    This class extends the IcontrolNetlist class and specializes in representing
    current-controlled current sources. It overrides the get_type_char method to return
    the current-controlled current source type character ('F').

    Parameters:
    - name (str): Name of the current-controlled current source netlist.
    - n1: Node 1 of the current-controlled current source netlist.
    - n2: Node 2 of the current-controlled current source netlist.
    - Vcontrol: I(Vcontrol) control the output current from the source 
    - gain: Gain associated with the current-controlled current source.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      current-controlled current source type character ('F').

    Note: The CurrentCurrentControlled inherits from the IcontrolNetlist class.
    """

    def __init__(self, name, n1, n2, Vcontrol, gain):
        """
        Constructor for a CurrentCurrentControlled instance with a name, nodes, a
        controlling current source, and a gain.
        :param name (str): Name of the current-controlled current source netlist.
        :param n1: Node 1 of the current-controlled current source netlist.
        :param n2: Node 2 of the current-controlled current source netlist.
        :param Vcontrol: I(Vcontrol) control the value of the current-controlled current source
        :param gain: Gain associated with the current-controlled current source.
        """
        super(CurrentCurrentControlled, self).__init__(name, n1, n2, Vcontrol, gain)

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The current-controlled current source type character ('F').
        """
        return 'F'

    
class MutualInductance(DipoleNetlist):
    """
    Represents a mutual inductance netlist for electric circuit simulation in spice format.

    This class extends the DipoleNetlist class and specializes in representing
    mutual inductances. It overrides the get_type_char method to return
    the mutual inductance type character ('K').

    Parameters:
    - name (str): Name of the mutual inductance netlist.
    - value: Value associated with the mutual inductance netlist.
    - L1: First inductance linked to L2 thanks to the mutual inductance netlist.
    - L2: Second inductance linked to L1 thanks to the mutual inductance netlist.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      mutual inductance type character ('K').

    Note: The MutualInductance inherits from the DipoleNetlist class.
    """

    def __init__(self, name, L1, L2, value):
        """
        Constructor for a MutualInductance instance with a name, value, and nodes.
        :param name (str): Name of the mutual inductance netlist.
        :param value: Value associated with the mutual inductance netlist.
        :param L1: First inductance linked to L2 thanks to the mutual inductance netlist.
        :param L2: Second inductance linked to L1 thanks to the mutual inductance netlist.
        """
        super(MutualInductance, self).__init__(name, L1, L2, value)

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The mutual inductance type character ('K').
        """
        return 'K'


class DummyNetlist(object):
    """
    Represents a dummy netlist for electric circuit simulation in spice format.

    This versatile class serves as a placeholder or dummy netlist, providing a simple
    character representation without actual circuit functionality. It can be employed
    for various purposes, such as adding commentary, representing complex expressions,
    or marking not-yet-implemented netlist elements.

    Parameters:
    - printed_char (str): Character to be printed as the dummy netlist.

    Methods:
    - _char_netlist(): Returns the character to be printed as the dummy netlist.
    """

    def __init__(self, printed_char):
        """
        Constructor for a DummyNetlist instance with a printed character.
        :param printed_char (str): Character to be printed as the dummy netlist.
        """
        self.printed_char = printed_char
    
    def _char_netlist(self):
        """
        Returns the character to be printed as the dummy netlist.
        """
        return self.printed_char


class SubCircuit(Netlist):
    """
    Represents a subcircuit netlist for electric circuit simulation in spice format.

    This class extends the Netlist class and specializes in representing
    subcircuits. It overrides the get_type_char method to return
    the subcircuit type character ('X').

    Parameters:
    - name (str): Name of the subcircuit netlist.
    - list_port (list): List of ports associated with the subcircuit.
    - subcircuit_name (str): Name of the subcircuit being referenced.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      subcircuit type character ('X').
    - _char_netlist(): Overrides the method in the base class. Returns a string
      representation of the subcircuit netlist in the spice format.

    Note: The SubCircuit inherits from the Netlist class.
    """

    def __init__(self, name, list_port, subcircuit_name):
        """
        Constructor for a SubCircuit instance with a name, list of ports, and
        the name of the referenced subcircuit.
        :param name (str): Name of the subcircuit netlist.
        :param list_port (list): List of ports associated with the subcircuit.
        :param subcircuit_name (str): Name of the subcircuit being referenced.
        """
        super(SubCircuit, self).__init__(name, subcircuit_name)
        self.list_port = list_port

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The subcircuit type character ('X').
        """
        return "X"

    def _char_netlist(self):
        """
        Overrides the method in the base class.
        :return: A string representation of the subcircuit netlist in the spice format.
        """
        net_char = self.full_name() + " " + " ".join(self.list_port) + " " + str(self.value)
        return net_char


class SubCircuitDefinition(Netlist):
    """
    Represents a subcircuit definition netlist for electric circuit simulation in spice format.

    This class extends the Netlist class and specializes in defining subcircuits.
    It overrides the get_type_char method to return the subcircuit definition type
    character ('.subckt').

    Parameters:
    - name (str): Name of the subcircuit definition netlist.
    - list_instance (list): List of instances associated with the subcircuit definition.
    - list_port (list): List of ports associated with the subcircuit definition.

    Methods:
    - get_type_char(): Overrides the method in the base class. Returns the
      subcircuit definition type character ('.subckt').
    - _char_netlist(): Overrides the method in the base class. Returns a string
      representation of the subcircuit definition netlist in the spice format.
    - generate_instance(name, list_port): Generates an instance of the subcircuit
      based on the defined format.

    Note: The SubCircuitDefinition inherits from the Netlist class.
    """

    def __init__(self, name, list_instance, list_port):
        """
        Constructor for a SubCircuitDefinition instance with a name, list of instances,
        and list of ports.
        :param name (str): Name of the subcircuit definition netlist.
        :param list_instance (list): List of instances to be placed within the subcircuit definition.
        :param list_port (list): List of ports associated with the subcircuit definition.
        """
        super(SubCircuitDefinition, self).__init__(name, name)  # Assuming subcircuit definition value is its name
        self.list_instance = list_instance
        self.list_port = list_port

    def get_type_char(self):
        """
        Overrides the method in the base class.
        :return: The subcircuit definition type character ('.subckt').
        """
        return ".subckt "

    def _char_netlist(self):
        """
        Overrides the method in the base class.
        :return: A string representation of the subcircuit definition netlist in the spice format.
        """
        net_char = self.full_name() + " " + " ".join(self.list_port) + "\n"
        for net in self.list_instance:
            net_char += net._char_netlist() + "\n"
        net_char += ".ends"
        return net_char

    def generate_instance(self, name, list_port):
        """
        Generates an instance of the subcircuit based on the defined format.
        :param name (str): Name of the generated subcircuit instance.
        :param list_port (list): List of ports associated with the generated subcircuit instance.
        :return: An instance of the SubCircuit class.
        :raises IOError: If the provided port list does not match the defined format.
        """
        if len(list_port) != len(self.list_port):
            raise IOError("Port list does not match with the defined format. "
                          +str(len(list_port)) + " given, " + str(len(self.list_port)) + " needed")
        return SubCircuit(name, list_port, self.name)



