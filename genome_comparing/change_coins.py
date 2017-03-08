#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


class MinCoins(object):
    """
    The MinCoins class uses objects to keep track of the dynamic
    programming for the mininum number of coins problem given a set of
    coins and an amount needed to make change.
    """

    def __init__(self, change, coins):
        self.change = int(change)
        self.coins = coins if isinstance(coins, dict) else MinCoins.coin_set_from_str(coins)
        self.min_coins = {0: 0}
        self.min_num = None

    @staticmethod
    def coin_set_from_str(coins_list):
        return dict((int(k), True) for k in re.split(',', coins_list))

    def get_min_coins(self):
        """
        Call the dynamic programming process to find the mininum number
        of coins for the change and coin set. Stores the result as
        min_coins for future recall without having to recalculate.
        :return: mininum coins for the object.
        """
        if not self.min_num:
            self.min_num = self.find_min_coins()
        return self.min_num

    def find_min_coins(self):
        """
        Dynamic programming function for calculating the minimum number
        of coins.
        :return: minimum coins for the current change
        """
        m = 1
        while m <= self.change:
            self.min_coins[m] = sys.maxsize
            for coin in self.coins:
                if coin <= m:
                    if self.min_coins[m - coin] + 1 < self.min_coins[m]:
                        self.min_coins[m] = self.min_coins[m - coin] + 1
            m += 1
        return self.min_coins[self.change]


# print(MinCoins(*sys.stdin.read().splitlines()).get_min_coins())