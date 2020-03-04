import numpy as np
from collections import defaultdict
from src_ours import *

class MCTS_Node(object):
    def __init__(self, state: planets_environment, parent = None):
        self._number_of_visits = 0.
        self._results = defaultdict(int)
        self.state = state
        self.parent = parent
        self.children = []







    def expand(self):
        action = self.untried_actions.pop()
        next_state = self.state.move(action)
        child_node = MCTS_Node(next_state, parent=self)
        self.children.append(child_node)
        return child_node

    def is_terminal_node(self):
        return self.state.is_game_over()