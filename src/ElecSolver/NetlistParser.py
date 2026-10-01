import ast
import operator
import re
import os
from numpy import sqrt, array, concatenate, arange
from .TemporalSystemBuilder import TemporalSystemBuilder

class NetlistParser():
    """Parse a SPICE-like netlist into component, node, and coupling mappings.

    Passive component values may be numeric SI literals or arithmetic expressions
    enclosed in braces. Parameter expressions are resolved in declaration order and
    may reference parameters declared earlier in the netlist.

    Attributes
    ----------
    file_path : str
        Path to the source netlist.
    node_map : dict[str, int]
        Mapping from netlist node names to integer node indices.
    dipole_map : dict[str, dict]
        Mapping from component names to mapped node indices. Passive entries also
        contain their resolved numeric value.
    resistors : dict[str, dict]
        Resistors keyed by name, with ``n1``, ``n2``, and raw ``value`` fields.
    inductors : dict[str, dict]
        Inductors keyed by name, with ``n1``, ``n2``, and raw ``value`` fields.
    capacitors : dict[str, dict]
        Capacitors keyed by name, with ``n1``, ``n2``, and raw ``value`` fields.
    couplings : dict[str, dict]
        Magnetic couplings keyed by name.
    real_couplings : dict[str, dict]
        Resistive couplings keyed by name.
    current : dict[str, dict]
        Independent and behavioral current sources keyed by name, with ``n1`` and
        ``n2`` fields.
    voltage : dict[str, dict]
        Independent and behavioral voltage sources keyed by name, with ``n1`` and
        ``n2`` fields.
    params : dict[str, float]
        Resolved parameter values keyed by parameter name.
    coupling_map : dict[str, dict]
        Magnetic coupling coordinates and resolved values keyed by coupling name.
    real_coupling_map : dict[str, dict]
        Resistive coupling coordinates and resolved values keyed by coupling name.
    max_index_node : int
        Next available node index after mapping.
    """
    # Regex patterns
    RESISTOR_PATTERN  = r'^([Rr][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    INDUCTOR_PATTERN  = r'^([Ll][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    CAPACITOR_PATTERN = r'^([Cc][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    COUPLING_PATTERN  = r'^([Kk][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    REALCOUPLING_PATTERN  = r'^([Ww][\w]*)\s+(\S+)\s+(\S+)\s+(.+)$'
    CURRENT_PATTERN = r'^([Ii][\w]*)\s+(\S+)\s+(\S+)(?:\s+.+)?$'
    VOLTAGE_PATTERN = r'^([Vv][\w]*)\s+(\S+)\s+(\S+)(?:\s+.+)?$'
    BEHAVIORAL_PATTERN = r'^([Bb][\w]*)\s+(\S+)\s+(\S+)\s+([IiVv])\s*='
    PARAM_PATTERN = r'^\.param\s+([\w]+)\s*=\s*(.+)$'
    SI_VALUE_PATTERN = r'\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)([a-zA-Zµ]*)\s*'
    SI_EXPRESSION_VALUE_PATTERN = (
        r'(?<![\w.])(?:\d+(?:\.\d*)?|\.\d+)'
        r'(?:[eE][+-]?\d+)?[a-zA-Zµ]+(?![\w])'
    )

    ARITHMETIC_OPERATORS = {
        ast.Add: operator.add,
        ast.Sub: operator.sub,
        ast.Mult: operator.mul,
        ast.Div: operator.truediv,
        ast.Pow: operator.pow,
        ast.Mod: operator.mod,
    }
    UNARY_OPERATORS = {
        ast.UAdd: operator.pos,
        ast.USub: operator.neg,
    }

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
        """Initialize a parser for a netlist file.

        Parameters
        ----------
        file_path : str
            Path to a ``.net``, ``.cir``, or ``.sp`` file.

        Raises
        ------
        FileNotFoundError
            If ``file_path`` does not exist.
        ValueError
            If the file extension is not supported.
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
        self.current = {}
        self.voltage = {}
        self.max_index_node = 0


    def map_netlist(self):
        """Parse the netlist and populate node, dipole, and coupling mappings."""
        self._parse_netlist()
        self._map_nodes()
        self._map_couplings()
        #print(self.dipole_map, self.coupling_map)


    def _parse_netlist(self):
        """Extract components, sources, couplings, and parameters from the netlist.

        Component values remain as raw strings until node mapping. Parameters are
        resolved immediately in declaration order.

        Raises
        ------
        ValueError
            If a parameter value or expression cannot be resolved.
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
        self.current = {m[0]: {'n1': m[1], 'n2': m[2]}
                for m in re.findall(self.CURRENT_PATTERN, data, re.MULTILINE)}
        self.voltage = {m[0]: {'n1': m[1], 'n2': m[2]}
                for m in re.findall(self.VOLTAGE_PATTERN, data, re.MULTILINE)}
        for name, n1, n2, source_type in re.findall(
                self.BEHAVIORAL_PATTERN, data, re.MULTILINE):
            sources = self.current if source_type.lower() == 'i' else self.voltage
            sources[name] = {'n1': n1, 'n2': n2}

        param_definitions = re.findall(self.PARAM_PATTERN, data, re.MULTILINE)
        self._parse_param_values(param_definitions)

    def _parse_param_values(self, param_definitions):
        """Resolve parameter definitions in declaration order.

        Parameters
        ----------
        param_definitions : iterable of tuple[str, str]
            Parameter names and raw expressions in netlist declaration order.

        Raises
        ------
        ValueError
            If an expression is invalid or references an unknown or later parameter.
        """
        self.params = {}
        for name, expression in param_definitions:
            try:
                self.params[name] = self._parse_si_value(expression)
            except ValueError as error:
                raise ValueError(
                    f"Invalid value for parameter '{name}': {error}"
                ) from error

    def _parse_si_literal(self, value_str):
        """Convert a numeric literal with an optional SI prefix to base units.

        Parameters
        ----------
        value_str : str
            Numeric literal, such as ``"1.5e-3"``, ``"10k"``, or ``"2u"``.

        Returns
        -------
        float
            Value converted to base units.

        Raises
        ------
        ValueError
            If the literal or SI prefix is invalid.
        """
        try:
            return float(value_str)
        except ValueError:
            pass

        match = re.fullmatch(self.SI_VALUE_PATTERN, value_str)
        if not match:
            raise ValueError(f"Invalid SI value format: '{value_str}'")

        number, prefix = match.groups()
        if prefix in self.SI_COEF:
            multiplier = self.SI_COEF[prefix]
        elif prefix.lower() in self.SI_COEF:
            multiplier = self.SI_COEF[prefix.lower()]
        else:
            raise ValueError(f"Unknown SI prefix: '{prefix}'")

        return float(number) * multiplier

    def _evaluate_param_expression(self, expression):
        """Evaluate arithmetic using already resolved parameters.

        Supported expressions contain real numeric literals, SI-prefixed literals,
        parameter names, parentheses, and the configured unary or binary arithmetic
        operators. Function calls and other Python syntax are rejected.

        Parameters
        ----------
        expression : str
            Expression without outer braces.

        Returns
        -------
        float
            Evaluated value in base units.

        Raises
        ------
        ValueError
            If the expression is empty, malformed, unsupported, non-numeric, or
            references an unknown or later parameter.
        """
        expression = expression.strip()
        if not expression:
            raise ValueError("Parameter expression cannot be empty")

        def replace_si_value(match):
            return repr(self._parse_si_literal(match.group(0)))

        normalized_expression = re.sub(
            self.SI_EXPRESSION_VALUE_PATTERN,
            replace_si_value,
            expression,
        )
        try:
            syntax_tree = ast.parse(normalized_expression, mode='eval')
        except SyntaxError as error:
            raise ValueError(
                f"Invalid parameter expression: '{expression}'"
            ) from error

        def evaluate(node):
            if isinstance(node, ast.Expression):
                return evaluate(node.body)
            if isinstance(node, ast.Constant):
                if isinstance(node.value, bool) or not isinstance(
                        node.value, (int, float)):
                    raise ValueError(
                        f"Unsupported value in parameter expression: '{expression}'"
                    )
                return node.value
            if isinstance(node, ast.Name):
                if node.id not in self.params:
                    raise ValueError(
                        f"Unknown or forward parameter reference '{node.id}'"
                    )
                return self.params[node.id]
            if isinstance(node, ast.BinOp):
                operation = self.ARITHMETIC_OPERATORS.get(type(node.op))
                if operation is None:
                    raise ValueError(
                        f"Unsupported operator in parameter expression: '{expression}'"
                    )
                return operation(evaluate(node.left), evaluate(node.right))
            if isinstance(node, ast.UnaryOp):
                operation = self.UNARY_OPERATORS.get(type(node.op))
                if operation is None:
                    raise ValueError(
                        f"Unsupported operator in parameter expression: '{expression}'"
                    )
                return operation(evaluate(node.operand))
            raise ValueError(
                f"Unsupported syntax in parameter expression: '{expression}'"
            )

        try:
            result = evaluate(syntax_tree)
        except (ZeroDivisionError, OverflowError) as error:
            raise ValueError(
                f"Cannot evaluate parameter expression '{expression}': {error}"
            ) from error

        if isinstance(result, bool) or not isinstance(result, (int, float)):
            raise ValueError(
                f"Parameter expression is not a real number: '{expression}'"
            )
        return float(result)

    def _parse_si_value(self, value_str):
        """Convert an SI literal or parameter expression to base units.

        Parameters
        ----------
        value_str : str
            Numeric SI literal, parameter name, or arithmetic expression. Expressions
            may optionally be enclosed in braces.

        Returns
        -------
        float
            Parsed or evaluated value in base units.

        Raises
        ------
        ValueError
            If the value is malformed, uses an unsupported SI prefix or expression
            syntax, or references an unknown parameter.
        """
        value_str = value_str.strip()
        if value_str.startswith("{") or value_str.endswith("}"):
            if not (value_str.startswith("{") and value_str.endswith("}")):
                raise ValueError(f"Unmatched expression braces: '{value_str}'")
            return self._evaluate_param_expression(value_str[1:-1])

        try:
            return self._parse_si_literal(value_str)
        except ValueError:
            pass
        return self._evaluate_param_expression(value_str)

    def _map_nodes(self):
        """Map component nodes and populate ``dipole_map``.

        Ground is always mapped to index zero. Passive entries contain mapped nodes
        and a resolved value; current and voltage source entries contain mapped nodes
        only.
        """
        # Create a mapping of node names to unique integers
        self.node_map = {'0': 0}  # Ground node
        self.dipole_map = {}
        passive_list = [self.resistors, self.inductors, self.capacitors]
        source_list = [self.current, self.voltage]
        dipole_list = passive_list + source_list
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
        for component in passive_list:
            for dipole in component:
                self.dipole_map[dipole] = {"nodes":[self.node_map[component[dipole]['n1']],
                                          self.node_map[component[dipole]['n2']]],
                                          "value":self._parse_si_value(component[dipole]['value'])}
        for component in source_list:
            for dipole in component:
                self.dipole_map[dipole] = {"nodes":[self.node_map[component[dipole]['n1']],
                                          self.node_map[component[dipole]['n2']]]}

    def _map_couplings(self):
        """Map valid magnetic and resistive couplings to inductor coordinates."""
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
        """Append passive component coordinates and values to solver arrays.

        Parameters
        ----------
        indexes : numpy.ndarray
            Node coordinate array with shape ``(2, n)``.
        values : numpy.ndarray
            Component value array with shape ``(n,)``.
        dipole_list : collections.abc.Mapping
            Mapping whose keys are passive component names in ``dipole_map``.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            Coordinate and value arrays with the requested components appended.
        """
        for dipole in dipole_list:
            nodes = self.dipole_map[dipole]["nodes"]
            value = self.dipole_map[dipole]["value"]
            indexes = concatenate((indexes, [[nodes[0]],[nodes[1]]]), axis=1)
            values = concatenate((values, [value]), axis=0)
        return indexes, values

    def _fill_array_coupling(self,indexes, values, coupling_map):
        """Append coupling coordinates and values to solver arrays.

        Parameters
        ----------
        indexes : numpy.ndarray
            Coupling coordinate array with shape ``(2, n)``.
        values : numpy.ndarray
            Coupling value array with shape ``(n,)``.
        coupling_map : collections.abc.Mapping
            Couplings keyed by name, each containing ``L_coords`` and ``value``.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            Coordinate and value arrays with the couplings appended.
        """
        for coupling in coupling_map:
            nodes = coupling_map[coupling]["L_coords"]
            value = coupling_map[coupling]["value"]
            indexes = concatenate((indexes, [[nodes[0]],[nodes[1]]]), axis=1)
            values = concatenate((values, [value]), axis=0)
        return indexes, values

    def generate_temporal_system(self):
        """Build a temporal system from the mapped passive components and couplings.

        ``map_netlist`` must be called before this method. Parsed current and voltage
        sources are metadata only and are not added to the generated system.

        Returns
        -------
        TemporalSystemBuilder
            Temporal system containing the parsed passive components and couplings.
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
