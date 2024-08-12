from ortools.sat.python import cp_model
import math
import sys
import time

serpent_dlct = [
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 0, -4, 0, -4, -4, 0, 4, 0, -4, 0, 0, 0, 4, 0, 0],
    [8, 0, 0, 0, -4, 0, 0, -4, -8, 0, 0, 0, 4, 0, 0, 4],
    [8, -4, 0, 0, 4, -4, 0, -4, 0, 0, -4, 0, 0, 4, 0, 0],
    [8, 0, 0, -8, 0, 0, 0, 0, -8, 0, 0, 8, 0, 0, 0, 0],
    [8, 4, 0, 0, 0, 0, -4, 0, 0, 0, 4, 0, -4, 0, -4, -4],
    [8, -4, -4, 0, 0, 0, 0, 0, 8, -4, -4, 0, 0, 0, 0, 0],
    [8, 0, 4, 0, 0, 0, -4, 0, 0, 4, 0, 0, -4, 0, -4, -4],
    [8, -4, 0, 0, -4, 0, -4, 4, 0, 0, -4, 0, 0, 0, 4, 0],
    [8, 0, 0, -8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 0, -4, 0, 4, 0, -4, -4, 0, -4, 0, 0, 0, 0, 4, 0],
    [8, 0, 0, 0, -4, 0, 0, -4, 0, 0, 0, -8, 4, 0, 0, 4],
    [8, 0, 4, 0, 0, -4, 0, 0, 0, 4, 0, 0, -4, -4, 0, -4],
    [8, -4, -4, 8, 0, 4, 4, 0, 0, -4, -4, 0, 0, -4, -4, 0],
    [8, 4, 0, 0, 0, -4, 0, 0, 0, 0, 4, 0, -4, -4, 0, -4],
    [8, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, -8, 0, -4, -4, 0]
];

chi_act = [[8, 8, 8, 8, 8, 8, 8, 8],
           [8, -8, 0, 0, 0, 0, 0, 0],
           [8, 0, -8, 0, 0, 0, 0, 0],
           [8, 0, 0, -8, 0, 0, 0, 0],
           [8, 0, 0, 0, -8, 0, 0, 0],
           [8, 0, 0, 0, 0, -8, 0, 0],
           [8, 0, 0, 0, 0, 0, -8, 0],
           [8, 0, 0, 0, 0, 0, 0, -8]];

serpent_act = [
    [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16],
    [16, 0, -8, 0, -8, -8, 0, 8, 0, -8, 0, 0, 0, 8, 0, 0],
    [16, 0, 0, 0, -8, 0, 0, -8, -16, 0, 0, 0, 8, 0, 0, 8],
    [16, -8, 0, 0, 8, -8, 0, -8, 0, 0, -8, 0, 0, 8, 0, 0],
    [16, 0, 0, -16, 0, 0, 0, 0, -16, 0, 0, 16, 0, 0, 0, 0],
    [16, 8, 0, 0, 0, 0, -8, 0, 0, 0, 8, 0, -8, 0, -8, -8],
    [16, -8, -8, 0, 0, 0, 0, 0, 16, -8, -8, 0, 0, 0, 0, 0],
    [16, 0, 8, 0, 0, 0, -8, 0, 0, 8, 0, 0, -8, 0, -8, -8],
    [16, -8, 0, 0, -8, 0, -8, 8, 0, 0, -8, 0, 0, 0, 8, 0],
    [16, 0, 0, -16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [16, 0, -8, 0, 8, 0, -8, -8, 0, -8, 0, 0, 0, 0, 8, 0],
    [16, 0, 0, 0, -8, 0, 0, -8, 0, 0, 0, -16, 8, 0, 0, 8],
    [16, 0, 8, 0, 0, -8, 0, 0, 0, 8, 0, 0, -8, -8, 0, -8],
    [16, -8, -8, 16, 0, 8, 8, 0, 0, -8, -8, 0, 0, -8, -8, 0],
    [16, 8, 0, 0, 0, -8, 0, 0, 0, 0, 8, 0, -8, -8, 0, -8],
    [16, 0, 0, 0, 0, 8, 8, 0, 0, 0, 0, -16, 0, -8, -8, 0],
];
inv3_act = [
    [8, 8, 8, 8, 8, 8, 8, 8],
    [8, 0, 0, 0, 0, 0, 0, -8],
    [8, 0, -8, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 0, -8, 0, 0],
    [8, -8, 0, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 0, 0, -8, 0],
    [8, 0, 0, -8, 0, 0, 0, 0],
    [8, 0, 0, 0, -8, 0, 0, 0]
];

inv5_act = [
    [32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32],
    [32, 0, -8, 0, 8, 0, 0, 0, -8, 8, 8, 0, 8, 0, -8, -8, 8, 0, 8, -8, 0, 0, 0, -8, -8, 0, 0, 0, -8, -8, 0, -8],
    [32, -8, 0, 0, -8, 8, 0, -8, 8, 8, 0, -8, -8, 0, -8, -8, 0, 0, 8, 0, 8, 0, 8, -8, 0, -8, 0, 0, 0, 0, -8, 0],
    [32, 0, -8, -8, 0, 0, 0, 0, -8, 0, -8, 0, -8, 0, 0, 8, -8, 0, 0, -8, -8, 0, 8, 0, 8, -8, 8, 8, -8, 8, 0, 0],
    [32, 0, 8, -8, 8, 0, 0, -8, 0, 8, 0, -8, 0, 0, 0, 0, -8, 0, -8, 0, 8, -8, -8, -8, 0, 0, 8, 8, -8, 0, 0, -8],
    [32, 0, -8, -8, 0, -8, 8, 0, 8, 0, -8, 0, -8, 8, -8, -8, 0, 0, -8, 8, 8, 0, 0, 8, -8, 0, -8, 0, 0, 0, 0, 0],
    [32, -8, 0, 0, -8, -8, 0, 8, -8, 0, 0, 0, 8, 8, 8, 0, 0, -8, 0, 0, 0, 0, -8, 0, 0, -8, -8, 8, -8, 8, -8, 0],
    [32, 8, -8, 0, -8, 0, 0, 8, 0, -8, 0, -8, 0, 8, 0, 8, 0, 0, -8, -8, 0, 0, -8, -8, 0, 0, 0, 0, 8, -8, -8, 8],
    [32, 8, 0, -8, 0, 0, 0, 0, -8, -8, -8, -8, 0, 8, 0, -8, 0, -8, 8, 0, 8, -8, 0, 0, 0, 0, 8, -8, 0, 8, -8, 0],
    [32, -8, 8, -8, 8, 0, 0, 0, 8, 0, -8, 8, 0, -8, 0, 0, 8, 8, 0, -8, 0, 0, -8, 0, 0, 0, 0, -8, -8, -8, -8, 0],
    [32, -8, -8, 0, 8, -8, 8, -8, 0, -8, 0, 8, -8, -8, 0, 0, 0, -8, 0, 8, 0, 0, -8, -8, 0, 8, 8, 0, 0, 0, 0, 0],
    [32, 8, 8, 0, 0, -8, 0, 8, 0, 0, 8, -8, -8, -8, 8, -8, 8, 0, 0, -8, -8, 0, -8, 0, -8, -8, 0, 0, 0, 0, 0, 0],
    [32, 0, -8, 8, -8, 0, 8, 0, 0, 0, 0, 0, 0, -8, 8, 0, -8, 0, -8, 0, 0, 0, 8, 8, -8, 0, 0, -8, -8, 8, -8, -8],
    [32, -8, 0, 0, 0, 0, 8, 0, -8, 8, 8, 0, -8, 0, 0, 0, -8, -8, 0, -8, 0, -8, 0, 0, 0, 8, -8, -8, 8, -8, 8, 0],
    [32, -8, 0, 8, 0, 0, 8, 8, 0, -8, 0, -8, 0, 0, -8, 8, 8, 0, -8, 0, -8, -8, 0, 0, 0, -8, 0, -8, 0, 0, 8, -8],
    [32, 0, 0, -8, 0, 0, 0, -8, -8, 0, 8, 8, 8, 0, 8, -8, 0, 0, -8, 0, -8, 8, 0, -8, 0, -8, -8, -8, 8, 0, 0, 0],
    [32, 0, 0, 0, -8, -8, 8, -8, 0, 8, -8, 0, 0, 8, 8, 0, 8, -8, 0, 0, -8, -8, 0, 0, -8, 0, 8, 0, 0, -8, 0, -8],
    [32, 0, 8, 0, -8, -8, -8, 0, -8, 0, 0, 0, -8, 0, -8, 8, 0, 0, -8, 0, 0, 0, 0, -8, 0, 8, 8, -8, -8, 0, 8, 8],
    [32, 8, 0, 0, 8, -8, -8, 0, 8, 0, 0, 0, 0, 0, -8, 0, -8, -8, 8, 0, 0, 8, 0, 0, 8, -8, 0, -8, 0, -8, -8, -8],
    [32, -8, 0, 0, 0, -8, -8, -8, 0, -8, 0, 0, 8, 0, 0, 0, -8, 8, 0, -8, 8, -8, 8, 0, -8, -8, 0, 8, 0, 0, 0, 8],
    [32, -8, -8, -8, 0, 0, -8, 0, 0, 0, 0, -8, 0, 8, 0, 0, -8, 0, 8, 8, -8, 8, -8, 0, -8, 8, 0, -8, 8, 0, 0, 0],
    [32, 0, 0, 0, -8, 0, 0, -8, 0, -8, 0, 8, 0, 0, -8, -8, -8, 0, 0, -8, -8, 8, -8, 8, 8, 8, 0, 0, 0, 8, 0, -8],
    [32, 8, -8, 8, 0, 8, -8, -8, 8, 0, 0, 0, -8, 0, 0, 0, 8, 0, 0, 0, 0, -8, -8, 8, 0, -8, -8, -8, -8, 0, 0, 0],
    [32, 8, 8, 0, 0, 8, -8, 0, 0, -8, -8, 0, -8, 0, 0, 8, 0, -8, 0, -8, 0, 8, 0, -8, -8, 0, -8, 0, 0, 8, 0, -8],
    [32, -8, 0, 0, 0, 0, -8, 0, -8, -8, 0, 8, -8, 0, 8, -8, 0, 8, -8, 8, 0, 0, 0, 8, 0, 0, 0, 8, 0, -8, -8, -8],
    [32, 0, 8, -8, 0, 8, -8, 0, 0, 0, 8, 8, 0, -8, -8, 0, 0, -8, -8, 0, -8, -8, 0, 0, -8, 0, 0, 8, 0, 0, -8, 8],
    [32, 0, 0, 0, -8, 8, 0, 0, -8, 0, -8, 0, 0, -8, -8, 0, -8, 0, 0, 8, 8, 0, -8, 0, -8, -8, 0, 0, 8, -8, 8, 8],
    [32, -8, -8, 8, 0, 0, -8, 0, 0, 8, -8, -8, 8, -8, 0, -8, 0, -8, -8, -8, 0, 0, 8, 0, 8, 0, 0, 0, 0, 0, 8, 0],
    [32, 0, 0, 8, 0, 0, 0, 8, 8, -8, -8, 0, 0, 0, 0, -8, -8, 8, 0, 8, -8, -8, 0, -8, 0, 0, -8, 0, -8, -8, 0, 8],
    [32, 0, 0, 8, 8, -8, 0, -8, -8, 0, -8, -8, 0, -8, -8, 8, 0, 8, 0, 0, -8, 0, 0, 0, 8, 8, -8, 0, 0, 0, -8, 0],
    [32, 0, 0, -8, -8, 8, 0, -8, 0, -8, 8, -8, 0, -8, 0, 0, 0, -8, 0, 0, 0, 8, 8, 8, 0, 0, -8, 0, -8, -8, 8, 0],
    [32, 0, -8, -8, -8, -8, -8, 8, 0, 0, 0, 0, 8, -8, 0, 0, 0, 8, 8, 0, 0, -8, 0, -8, 8, 0, -8, 0, 8, 0, 0, -8],
];

inv4_act = [[16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16],
            [16, -8, 0, -8, 0, -8, 8, 0, 0, 0, 0, 0, 0, 0, 8, -8],
            [16, -8, 0, 0, 0, 0, 0, -8, -8, 0, -8, 8, 0, 0, 0, 8],
            [16, 8, -8, 0, 0, 0, 0, 0, -8, -8, 0, 0, 0, 8, 0, -8],
            [16, 0, 0, -8, -8, 8, 0, 8, -8, 0, 0, 0, 0, -8, 0, 0],
            [16, 0, 0, 8, 0, -8, 0, 0, -8, 8, 0, -8, -8, 0, 0, 0],
            [16, 0, 0, 0, -8, 0, 0, -8, 8, -8, 0, 0, -8, 0, 8, 0],
            [16, 0, -8, 0, 0, 0, 8, 0, 8, 0, -8, -8, 0, -8, 0, 0],
            [16, -8, -8, 8, -8, 0, 0, 0, 0, 0, 8, 0, 0, 0, -8, 0],
            [16, 0, 0, -8, 0, 0, 0, -8, 0, 0, 0, -8, 8, 8, -8, 0],
            [16, 8, 0, 0, -8, -8, -8, 0, 0, 0, -8, 0, 8, 0, 0, 0],
            [16, -8, 8, 0, 0, 0, -8, 8, 0, -8, 0, -8, 0, 0, 0, 0],
            [16, 0, -8, -8, 8, 0, -8, 0, 0, 0, 0, 0, -8, 0, 0, 8],
            [16, 0, 0, 0, 0, 0, -8, -8, 0, 8, 8, 0, 0, -8, 0, -8],
            [16, 0, 0, 0, 8, -8, 0, 0, 0, -8, 0, 8, 0, -8, -8, 0],
            [16, 0, 8, 0, 0, 8, 0, 0, 0, 0, -8, 0, -8, 0, -8, -8]];

ascon_act = [
    [32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32],
    [32, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0],
    [32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 32, -32, -32, -32, -32, 32, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0],
    [32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, -32, 0, -32, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 32, 0, 32, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, -32, 0, -32, 0, -32],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, -32, 32, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 32, 0, 0, 0, 0, 0, -32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 32, 0, 0, 0, -32, 0, 0, 0, -32, 0, 0, 0, 0, -32, 0, 0, 0, -32, 0, 0, 0, 32, 0, 0, 0, 32, 0, 0],
    [32, -32, 0, 0, 0, 0, 0, 0, 32, -32, 0, 0, 0, 0, 0, 0, -32, 32, 0, 0, 0, 0, 0, 0, -32, 32, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 32, 0, 0, 0, 0, -32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 32, 32, 0, 0, 0, 0, -32, -32, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 32, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 32, 0, 0, 0, 0, -32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, -32, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, -32, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, -32, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, -32, 0, 32, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];

keccak_act = [
    [32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     32, 32, 32],
    [32, -32, 32, -32, 32, -32, 32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, -32, 0, 32, 0, -32, 0, 32, 0, -32, 0, 32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, -32, 0, 32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, -32, 0, 0, 0, 32, 0, 0, 0, -32, 0, 0, 0, 32, 0, 0, 0, -32, 0, 0, 0, 32, 0, 0, 0, -32, 0, 0, 0],
    [32, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, -32, 0, 0, 0, 32, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 32, 0, 0, 0, 0, 0, 0, -32, -32, 0, 0, 0, 0, 0, 0, 32, 32, 0, 0, 0, 0, 0, 0, -32, -32, 0, 0, 0, 0, 0, 0],
    [32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 32, 0, 0, 0, 0, -32, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 32, 32, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, -32, -32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, -32, 32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, -32, 0, 0, 0, 0, 0],
    [32, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, -32, 0],
    [32, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, -32, 0, 0, 0, 0],
    [32, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32, 0, 0, 0, 0, 0, 0, 0],
    [32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -32]];


gift_act = [[16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16],
            [16, 0, 0, 0, 0, 0, 0, 0, -8, 0, 0, -8, -8, 0, 0, 8],
            [16, 0, 0, -16, -8, 0, 0, 8, 0, 0, 0, 0, -8, 0, 0, 8],
            [16, 0, 0, 0, -8, 0, 0, 8, -8, 0, 0, -8, 0, 0, 0, 0],
            [16, -16, 0, 0, -8, 8, 0, 0, 8, -8, -8, 8, -8, 8, 0, 0],
            [16, 0, 0, 0, -8, 8, 0, 0, -8, 0, 0, -8, 8, 0, -8, 0],
            [16, 0, -16, 0, 8, 0, -8, 0, 8, -8, -8, 8, 8, 0, -8, 0],
            [16, 0, 0, 0, 8, 0, -8, 0, -8, 0, 0, -8, -8, 8, 0, 0],
            [16, -16, -16, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [16, 0, 0, 0, 0, 0, 0, 0, 0, -8, -8, 0, 0, -8, 8, 0],
            [16, 0, 0, -16, 0, -8, 8, 0, 0, 0, 0, 0, 0, -8, 8, 0],
            [16, 0, 0, 0, 0, -8, 8, 0, 0, -8, -8, 0, 0, 0, 0, 0],
            [16, 16, 0, 0, 0, 0, -8, -8, 0, 0, 0, 0, 0, 0, -8, -8],
            [16, 0, 0, 0, 0, 0, -8, -8, 0, 8, 8, 0, 0, -8, 0, -8],
            [16, 0, 16, 0, 0, -8, 0, -8, 0, 0, 0, 0, 0, -8, 0, -8],
            [16, 0, 0, 0, 0, -8, 0, -8, 0, 8, 8, 0, 0, 0, -8, -8]];

present_act = [[16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16],
               [16, -16, 0, 0, 0, 0, 0, 0, 0, 0, -16, 16, 0, 0, 0, 0],
               [16, 0, 0, -8, -8, 0, -8, 8, 0, -8, 0, 0, 0, 0, 0, 8],
               [16, 0, -8, 0, 0, -8, 0, 0, 8, 0, 0, 0, -8, -8, 8, 0],
               [16, 0, 0, -8, -8, 0, 0, 0, 0, -8, 0, 0, -8, 8, 0, 8],
               [16, 0, 8, 0, 0, -8, -8, -8, -8, 0, 0, 0, 0, 0, 8, 0],
               [16, 0, -8, 8, 0, 0, 0, 0, -8, 8, 0, -16, 0, 0, 0, 0],
               [16, 0, 0, 0, 0, 0, 8, -8, 0, 0, 0, -16, 8, -8, 0, 0],
               [16, -16, -8, 8, 0, 0, 0, 0, -8, 8, 0, 0, 0, 0, 0, 0],
               [16, 16, 0, 0, -8, -8, 0, 0, 0, 0, 0, 0, 0, 0, -8, -8],
               [16, 0, 0, -8, 0, 8, -8, 8, 0, -8, 0, 0, 0, 0, -8, 0],
               [16, 0, 8, 0, 8, 0, 0, 0, -8, 0, 0, 0, -8, -8, 0, -8],
               [16, 0, 0, -8, 0, 8, 0, 0, 0, -8, 0, 0, -8, 8, -8, 0],
               [16, 0, -8, 0, 8, 0, -8, -8, 8, 0, 0, 0, 0, 0, 0, -8],
               [16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -16, 0, 0, 0, 0, 0],
               [16, 0, 0, 0, -8, -8, 8, -8, 0, 0, 16, 0, 8, -8, -8, -8]];

skinny_act = [[16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16],
              [16, 0, 0, 0, 16, 0, 0, 0, -16, 0, 0, 0, -16, 0, 0, 0],
              [16, -8, 0, -8, 0, -8, 0, 8, 16, -8, 0, -8, 0, -8, 0, 8],
              [16, 0, 0, 0, 0, 0, 0, 0, -16, 0, 0, 0, 0, 0, 0, 0],
              [16, 0, -8, 0, 0, 0, -8, 0, 0, 8, -8, -8, 0, 8, 8, -8],
              [16, 8, -8, -8, 0, 8, -8, -8, 0, 0, -8, 0, 0, 0, 8, 0],
              [16, 0, 0, 0, 0, 0, 0, 0, 0, -8, 0, 8, 0, -8, 0, -8],
              [16, -8, 0, 8, 0, -8, 0, -8, 0, 0, 0, 0, 0, 0, 0, 0],
              [16, 0, 8, 0, -16, 0, -8, 0, 0, 0, 8, 0, 0, 0, -8, 0],
              [16, 0, 8, 0, -16, 0, -8, 0, 0, 0, 8, 0, 0, 0, -8, 0],
              [16, 0, 0, -8, 0, 0, 0, 8, 0, 0, 0, -8, -16, 0, 0, 8],
              [16, -8, 0, 0, 0, -8, 0, 0, 0, -8, 0, 0, 16, -8, 0, 0],
              [16, 8, -8, 0, 0, 8, 8, 0, 0, 0, -8, -8, 0, 0, -8, -8],
              [16, 0, -8, -8, 0, 0, 8, -8, 0, 8, -8, 0, 0, 8, -8, 0],
              [16, -8, 0, 0, 0, -8, 0, 0, 0, 0, 0, 8, 0, 0, 0, -8],
              [16, 0, 0, 8, 0, 0, 0, -8, 0, -8, 0, 0, 0, -8, 0, 0]];

class dlct_act_to_sbox():
    def __init__(self, n, workers, T):
        self.n = n
        self.workers = workers
        self.T = T
        self.model = cp_model.CpModel()
        self.solver = cp_model.CpSolver()

    def solutions_pairs(self):
        INDEXS = []
        for alpha in range(1, self.n):
            T = []
            for i in range(self.n):
                for j in range(i + 1, self.n):
                    if i ^ j == alpha:
                        T.append([i, j])
            INDEXS.append(T)
        return INDEXS

    def generate_model(self):
        N = int(math.log(self.n, 2))
        INDEXS = self.solutions_pairs()
        Y = [self.model.NewIntVar(0, self.n - 1, f'Y[{i}]') for i in range(self.n)]
        A = {}

        for i in range(self.n):
            for j in range(self.n):
                A[i, j] = self.model.NewBoolVar(f'A[{i},{j}]')

        X = {}
        for i in range(self.n):
            for j in range(self.n):
                X[i, j] = self.model.NewBoolVar(f'X[{i},{j}]')

        for j in range(self.n):
            self.model.Add(sum(X[i, j] for i in range(self.n)) == 1)

        for i in range(self.n):
            self.model.Add(sum(X[i, j] for j in range(self.n)) == 1)

        for i in range(self.n):
            self.model.Add(sum(j * X[i, j] for j in range(self.n)) == Y[i])

        for i in range(self.n):
            self.model.Add(Y[i] == sum((1 << j) * A[i, j] for j in range(N)))

        S = {}
        V = {}
        M = {}
        for row in range(self.n - 1):
            for la in range(self.n - 1):
                for i in range(int(self.n / 2)):
                    S[row, la, i] = self.model.NewBoolVar(f'S[{row},{la},{i}]')
                    M[row, la, i] = self.model.NewBoolVar(f'M[{row},{la},{i}]')
                    V[row, la, i] = self.model.NewIntVar(0, N, f'V[{row},{la},{i}]')

        T1 = [0] * N
        T2 = [0] * N

        for row in range(self.n-1):
            for la in range(self.n - 1):
                for i in range(int(self.n / 2)):
                    for k in range(N):
                        T1[k] = A[INDEXS[row][i][0], k]
                        T2[k] = A[INDEXS[row][i][1], k]
                    self.model.Add(M[row, la, i] == 1 - S[row, la, i])
                    temp = 0
                    for k in range(N):
                        temp += (((la + 1) >> k) & 1) * T1[k]
                        temp += (((la + 1) >> k) & 1) * T2[k]
                    temp += M[row, la, i]
                    self.model.Add(temp == 2 * V[row, la, i])

                self.model.Add(sum(S[row, la, j] for j in range(int(self.n / 2))) == int(
                    (self.T[row + 1][la + 1] / 2 + self.n / 2) / 2))            # In case of ACT
                #self.model.Add(sum(S[row, la, j] for j in range(int(self.n / 2))) == int((self.T[row + 1][la + 1]  + self.n / 2) / 2))  # In case of DLCT

        return Y

    def solve_model(self):
        Y = self.generate_model()
        self.solver.parameters.num_search_workers = self.workers
        # self.solver.parameters.max_time_in_seconds = 10.0
        # self.solver.parameters.log_search_progress = True

        status = self.solver.Solve(self.model)

        if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
            Z = [self.solver.Value(Y[i]) for i in range(self.n)]
            print(f'S = SBox({Z})')
            print(f'Time = {self.solver.WallTime()} s')
        else:
            print('No solution found.')
            print(f'Time = {self.solver.WallTime()} s')


if __name__ == '__main__':
    n = len(sys.argv)
    # 1: n
    # 2: threads
    D = dlct_act_to_sbox(int(sys.argv[1]), int(sys.argv[2]), inv5_act)
    D.solve_model()