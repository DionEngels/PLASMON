import _thread


def to_print():
    global x
    print(x)
    x *= 2
    print(x)


if __name__ == '__main__':
    x = 5

    _thread.start_new_thread(to_print, ())
