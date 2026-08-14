
from pathlib import Path
from matplotlib import pyplot as plt

def sliding_window_plot(values, path='test_plots/', title='plot', save=True, show=False):
    plt.figure()
    plt.plot(range(len(values)), values)
    plt.xlabel('Window index')
    plt.ylabel('Values')
    plt.title(title)
    plt.grid(True)
    plt.ylim(0, 100)
    plt.tight_layout()

    Path(path).mkdir(parents=True, exist_ok=True)
    if save:
        plt.savefig(f'{path}/{title}.png')
    if show:
        plt.show()

    plt.close()