import numpy as np
import time


def partition(arr, low, high):
    i = (low - 1)  # index of smaller element
    pivot = arr[high]  # pivot
    for j in range(low, high):
        if arr[j] <= pivot:
            i = i + 1
            arr[i], arr[j] = arr[j], arr[i]

    arr[i + 1], arr[high] = arr[high], arr[i + 1]
    return i + 1


def find_k_best(arr, low, high, k):
    if low < high:
        pi = partition(arr, low, high)
        if pi < k:
            return find_k_best(arr, pi, high, k)
        elif pi > k:
            return find_k_best(arr, low, pi - 1, k)
        else:
            return arr[pi-1]


if __name__ == '__main__':
    arr = np.arange(1, 100)
    np.random.shuffle(arr)
    k = 1
    print(find_k_best(arr, 0, len(arr) - 1, k))