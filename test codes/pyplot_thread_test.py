import matplotlib.pyplot as plt
import _thread


def not_thread_safe(x, y):
    plt.plot(x, y)
    plt.xlabel('x label')
    plt.ylabel('y label')
    plt.title("I don't get along with multi threading")
    plt.show()


def thread_safe(x, y):
    try:
        fig, ax = plt.subplots()  # fig : figure object, ax : Axes object
        ax.plot(x, y)
        ax.set_xlabel('x label')
        ax.set_ylabel('y label')
        ax.set_title("I am compatible with multithreading")
        plt.show()
    except Exception:
        import traceback
        print(traceback.format_exc())


if __name__ == '__main__':
    x = [1, 2, 3, 4, 5, 6]
    y = [7, 5, 9, 1, 4, 3]

    not_thread_safe(x, y)

    _thread.start_new_thread(thread_safe, (x, y))
