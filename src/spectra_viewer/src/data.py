from src.utils import create_peaks
from random import randint

initial_positive = create_peaks(randint(30, 60))
initial_negative = create_peaks(randint(30, 60))

rand_pep = randint(200, 1200)+randint(100, 999)/1000
rand_pepB = randint(201, 1199)+randint(100, 999)/1000
