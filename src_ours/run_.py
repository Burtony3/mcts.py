import numpy as np
from nodes import *
from search import MonteCarloTreeSearch
from tictactoe import TicTacToeGameState


def init():
    state = np.zeros((3, 3))