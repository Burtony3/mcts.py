import numpy as np
from nodes2 import *
from search import MonteCarloTreeSearch
from tictactoe2 import FlyByState
import spiceypy as spice


def init():
    state = np.array([3])
    initial_state = FlyByState(state=state, next_to_move=1, l=45)
    root = MonteCarloTreeSearchNode(state=initial_state, parent=None)
    mcts = MonteCarloTreeSearch(root)
    best_node = mcts.best_action(1)
    c_state = best_node.state
    c_board = c_state.board
    return c_state,c_board

#
# def graphics(board):
#     for i in range(3):
#         print("")
#         print("{0:3}".format(i).center(8)+"|", end='')
#         for j in range(3):
#             if c_board[i][j] == 0:
#                 print('_'.center(8), end='')
#             if c_board[i][j] == 1:
#                 print('X'.center(8), end='')
#             if c_board[i][j] == -1:
#                 print('O'.center(8), end='')
#     print("")
#     print("______________________________")
#
#
# def get_action(state):
#     try:
#         location = input("Your move: ")
#         if isinstance(location, str):
#             location = [int(n, 10) for n in location.split(",")]
#         if len(location) != 2:
#             return -1
#         x = location[0]
#         y = location[1]
#         move = TicTacToeMove(x, y, -1)
#     except Exception as e:
#         move = -1
#     if move == -1 or not state.is_move_legal(move):
#         print("invalid move")
#         move = get_action(state)
#     return move
#
#
# def judge(state):
#     if state.is_game_over():
#         if state.game_result == 1.0:
#             print("You lose!")
#         if state.game_result == 0.0:
#             print("Tie!")
#         if state.game_result == -1.0:
#             print("You Win!")
#         return 1
#     else:
#         return -1


c_state,c_board = init()
# graphics(c_board)
# #
# #
# while True:
#     move1 = get_action(c_state)
#     c_state = c_state.move(move1)
#     c_board = c_state.board
#     graphics(c_board)
#
#     board_state = FlyByState(state=c_board, next_to_move=1)
#     root = MonteCarloTreeSearchNode(state=board_state, parent=None)
#     mcts = MonteCarloTreeSearch(root)
#     best_node = mcts.best_action(1000)
#     c_state = best_node.state
#     c_board = c_state.board
#     graphics(c_board)
#     if judge(c_state)==1:
#         break
#     elif judge(c_state)==-1:
#         continue
