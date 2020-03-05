import numpy as np
from node import *
from search import MonteCarloTreeSearch
from sequenceEnv import sequenceState


def init():
    # Defines the Root of the Tree (Pre Launch @ Earth)
    state = np.array([3])
    # User inputs? I think this is where they should go or where we should define them somehow to be inputs for the
    # environment
    initial_state = sequenceState(state=state, next_to_move=1, l=45)
    root = MonteCarloTreeSearchNode(state=initial_state, parent=None)
    mcts = MonteCarloTreeSearch(root)
    best_node = mcts.best_action(1)
    c_state = best_node.state
    c_board = c_state.board
    return c_state, c_board


c_state, c_board = init()
