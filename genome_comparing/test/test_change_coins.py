from genome_comparing.change_coins import *


__author__ = 'Michael Lockwood'
__github__ = 'mlockwood'
__email__ = 'lockwm@uw.edu'


def test_coin_set_from_str():
    assert MinCoins.coin_set_from_str('50,25,20,10,5,1') == {50: True, 25: True, 20: True, 10: True, 5: True, 1: True}


def test_min_coins_class():
    assert MinCoins(40, MinCoins.coin_set_from_str('50,25,20,10,5,1')).get_min_coins() == 2