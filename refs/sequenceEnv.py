import numpy as np
import random


class mctsMove(object):
    def __init__(self, planet, value):
        self.planet = planet
        self.value = value


class sequenceState(object):
    def __init__(self, state, next_to_move=1, l=45):
        #Sequence
        self.sequence = state
        #Player to Move -> Planet or Epoch
        self.next_to_move = next_to_move
        #Grid Resolution Parameter
        self.l = l

    @property
    def game_result(self):
        # delta v function
        return None

    def is_game_over(self):
        # constraints
        # if arrival planet or if delta V cumulative cntraint is violated

        return self.game_result is not None

    def move(self, move):
        newSequence = np.copy(self.sequence)
        newSequence = np.append(newSequence, move.planet)
        return sequenceState(newSequence, -self.next_to_move)

    def get_legal_actions(self):
        if self.next_to_move > 0:
            # Periods of Planets (seconds)
            T = np.array([88.0, 224.7, 365.2, 687.0, 4331, 10747, 30589, 59800]) * 84600
            if 0 <= 1 < len(self.sequence):
                epochs = [self.sequence[1] + x * (self.l * int(T[int(self.sequence[-1])])) / 360 for x in range(9)]
                return [mctsMove(epoch, self.next_to_move) for epoch in epochs]
            else:
                epochs = [x * (self.l * int(T[self.sequence[-1]])) / 360 for x in range(9)]
                return [mctsMove(epoch, self.next_to_move) for epoch in epochs]
        else:
            # Possible Planet Selection
            planets = [1, 2, 3, 4, 5, 6, 7, 8]
            return [mctsMove(planet, self.next_to_move) for planet in planets]
