import math

def zaokrnapolovinu(broj):
    broj = math.floor(broj/0.5)
    broj = broj * 0.5
    return(broj)

def radianiustepene(broj):
    return(broj * 180 / math.pi)

def stepeniuradiane(broj):
    return(broj*math.pi/180)
