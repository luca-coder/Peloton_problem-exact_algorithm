import math
import timeit

import Algorithm
peloton_speed = 2
ralph_speed = 1

if __name__ == '__main__':

    start_time = timeit.default_timer()
    Algorithm.solve()
    end_time = timeit.default_timer()

    total_time = end_time - start_time
    print("total time")
    print(total_time)
