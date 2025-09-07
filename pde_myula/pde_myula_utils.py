import numpy as np
import cuqi

def sample_in_batches(sampler, Ns, Nt, A):

    N_batches = 10
    all_samples_list = []    
    for batch_idx in range(N_batches):
        print(f'batch number {batch_idx} of {N_batches}')
        _ = sampler.sample(int(Ns/N_batches))
        samples_batch = sampler.get_samples()
        all_samples_list.append(samples_batch.burnthin(Nb=0, Nt=Nt).samples)
        sampler._samples = []

    all_samples_array = np.concatenate(all_samples_list, axis=1)
    posterior_samples = cuqi.samples.Samples(
        all_samples_array, geometry=A.domain_geometry)

    return posterior_samples

def plot_figure_12(kappa_true, y_true, y_obs, posterior_samples, line_samples, exact_line, xx):

    import os
    from matplotlib import ticker
    import matplotlib.pyplot as plt
    
    # Set up matplotlib
    SMALL_SIZE = 7
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 9
    plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    
    # Use latex package
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{bm}')

    # Create the figure
    cm_to_in = 1/2.54
    fig, axs = plt.subplots(nrows=2, ncols=3,
                            figsize=(17.8*cm_to_in, 9.8*cm_to_in),
                            layout="constrained")

    # (a)
    plt.sca(axs[0,0])

    im = kappa_true.plot(subplots=False, vmin=-0.55, vmax=0.02, mode='color')
    inset_axes = plt.gca().inset_axes([1.04, 0.2, 0.05, 0.6])
    fig.colorbar(im[0], ax=plt.gca(), cax=inset_axes)
    plt.gca().set_ylim(0, 1)
    plt.gca().set_xlim(0, 1)

    plt.gca().set_title('(a) Exact solution')
    plt.ylabel('$x_2$')
    plt.gca().yaxis.labelpad = -5
    plt.xlabel('$x_1$')
    plt.gca().xaxis.labelpad =  2
    
    # (b)
    plt.sca(axs[0,1])
    im = y_true.plot(subplots=False)
    
    inset_axes = plt.gca().inset_axes([1.04, 0.2, 0.05, 0.6])
    fig.colorbar(im[0], ax=plt.gca(), cax=inset_axes)
    plt.ylabel('$x_2$')
    plt.gca().yaxis.labelpad = -5
    plt.gca().set_title('(b) Exact data')
    # control number of ticks of the colorbar
    inset_axes.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))

    plt.xlabel('$x_1$')
    plt.gca().xaxis.labelpad =  2# -5 -5
    
    
    # (c)
    plt.sca(axs[0,2])

    im = y_obs.plot(subplots=False)
    inset_axes = plt.gca().inset_axes([1.04, 0.2, 0.05, 0.6])
    fig.colorbar(im[0], ax=plt.gca(), cax=inset_axes)
    plt.ylabel('$x_2$')
    plt.gca().yaxis.labelpad = -5
    # control number of ticks of the colorbar
    inset_axes.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
    plt.xlabel('$x_1$')
    plt.gca().xaxis.labelpad =  2# -5 -5
    plt.gca().set_title('(c) Noisy data')
    
    # (d)
    plt.sca(axs[1,0])
    im = posterior_samples.plot_mean(
        subplots=False, vmin=-0.55, vmax=0.02, mode='color')
    inset_axes = plt.gca().inset_axes([1.04, 0.2, 0.05, 0.6])
    fig.colorbar(im[0], ax=plt.gca(), cax=inset_axes)
    plt.gca().set_ylim(0, 1)
    plt.ylabel('$x_1$')
    plt.gca().yaxis.labelpad = -5
    plt.gca().set_xlim(0, 1)
    
    plt.xlabel('$x_1$')
    plt.gca().xaxis.labelpad =  2# -5 -5
    plt.gca().set_title('(d) Posterior mean')
    
    # (e)
    plt.sca(axs[1,1])
    im = posterior_samples.funvals.vector.plot_std(subplots=False)

    inset_axes = plt.gca().inset_axes([1.04, 0.2, 0.05, 0.6])
    cb = fig.colorbar(im[0], ax=plt.gca(), cax=inset_axes)
    cb.locator = ticker.MaxNLocator(nbins=4)
    plt.ylabel('$x_2$')
    plt.gca().yaxis.labelpad = -5
    plt.xlabel('$x_1$')
    plt.gca().xaxis.labelpad = 2# -5
    plt.gca().set_title('(e) Posterior STD')
    
    # (f)
    plt.sca(axs[1,2])
    
    line_samples_obj = cuqi.samples.Samples(line_samples.T, geometry=cuqi.geometry.Continuous1D(xx))
    lines = line_samples_obj.plot_ci( exact=exact_line)
    plt.legend(ncol=2, frameon=False, loc=2)
    plt.xlim([0,1])
    plt.ylim([-1, 0.5])
    plt.xlabel('$x_1$')
    plt.gca().xaxis.labelpad = 2# -5
    plt.gca().set_title('(f) Posterior CI')
    plt.gca().set_box_aspect(1)


def get_samples_at_line(posterior_samples, xx, yy):
    line_samples = []

    samples_funvals = posterior_samples.funvals

    for i, fun_i in enumerate(samples_funvals):
        temp_list = []
        for xx_i, yy_i in zip(xx, yy):
            temp_list.append(fun_i(xx_i, yy_i))
        line_samples.append(temp_list)

    line_samples = np.array(line_samples)
    return line_samples