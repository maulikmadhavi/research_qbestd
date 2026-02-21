"""This module contains functions for plotting."""

import matplotlib.pyplot as plt

def plot_wave(t, wave, freq):
    plt.figure()
    plt.plot(t, wave)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(f'Sine Wave at {freq} Hz')
    plt.show()
    
    
if __name__ == '__main__':
    from gen_data import get_20hz_sine_wave, get_10hz_sine_wave
    x1, t1 = get_20hz_sine_wave()
    x2, t2 = get_10hz_sine_wave()
    f1, f2 = 20, 10
    plot_wave(t1, x1, f1)
    plot_wave(t2, x2, f2)
    