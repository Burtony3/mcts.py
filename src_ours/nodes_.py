import numpy as np
from collections import defaultdict
from src_ours import *


class MCTS_Node(object):
    def __init__(self, state: planets_environment, parent=None):
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

    def rollout(self):
        current_rollout_state = self.state
        while not current_rollout_state.is_game_over():
            possible_moves = current_rollout_state.get_legal_actions()
            action = self.rollout_policy(possible_moves)
            current_rollout_state = current_rollout_state.move(action)
        return current_rollout_state.game_result

    def rollout_policy(self, possible_moves):
        return possible_moves[np.random.randint(len(possible_moves))]