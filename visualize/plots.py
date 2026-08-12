
from pathlib import Path
from matplotlib import pyplot as plt

def plot_gc_profile(gc_values, sequence_name="", save=True, show=False):

    plt.figure()
    plt.plot(range(len(gc_values)), gc_values)
    plt.xlabel('Window index')
    plt.ylabel('GC values')
    plt.title(sequence_name)
    plt.grid(True)
    plt.ylim(0, 100)      # GC% всегда 0-100
    plt.tight_layout()

    Path("test_plots/").mkdir(parents=True, exist_ok=True)
    if save:
        plt.savefig('test_plots/gc_profile.png')
    if show:
        pass
    plt.close()