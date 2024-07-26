'''
A selection of rate functions, i.e., speed curves for animations.
(1) https://docs.manim.community/en/stable/reference/manim.utils.rate_functions.html
(2) https://easings.net/
'''
# Indicate which objects should be imported with "from spacing import *"
__all__ = ["spacing"]

import numpy as np

# =============================================================================
# Define Spacing Functions
# =============================================================================
def easeOutExpo(x, base=2, power=-10):
    return 1 - base**(power * x)

def easeOutPoly(x, n=1.5):
    return 1 - (1 - x)**n

def sigmoid(x: float) -> float:
    r"""Returns the output of the logistic function.

    The logistic function, a common example of a sigmoid function, is defined
    as :math:`\frac{1}{1 + e^{-x}}`.

    References
    ----------
    - https://en.wikipedia.org/wiki/Sigmoid_function
    - https://en.wikipedia.org/wiki/Logistic_function
    """
    return 1.0 / (1 + np.exp(-x))

def smooth(t: float, inflection: float = 10.0) -> float:
    error = sigmoid(-inflection / 2)
    err = (sigmoid(inflection * (t - 0.5)) - error) / (1 - 2 * error)
    return np.minimum(np.maximum(err, 0), 1)

funcs = {
        'circle'        : lambda x, n=2: (1 - (x - 1)**n)**(1/n),
        'easeOutExpo'   : easeOutExpo,
        'poly'          : easeOutPoly,
        'root'          : lambda x, n=2: x**(1/n),
        'linear'        : lambda x: easeOutPoly(x, n=1),
        'quadratic'     : lambda x: easeOutPoly(x, n=2),
        'cubic'         : lambda x: easeOutPoly(x, n=3),
        # 'sigmoid'       : smooth,
    }


def spacing(mi, ma, num, func_name='linear', reverse=False, log_scale=False, **kwargs):
    if log_scale:
        # mi, ma = np.log10(mi), np.log10(ma)
        mi, ma = 10**mi, 10**ma
    x = np.linspace(0, 1, num)
    y = evaluate_function(x, func_name, **kwargs)
    y = (1 - y)[::-1] if reverse else y
    # result = mi + (ma - mi) * y
    # return 10**result if log_scale else result
    return mi + (ma - mi) * y



# =============================================================================
# MAIN
# =============================================================================
def main():
    trim = False

    # Choose what to plot
    plot_bar_chart(trim=trim)
    # plot_specific_function('all', trim=trim)
    pass


# =============================================================================
# OPTIONS
# =============================================================================
def plot_bar_chart(trim=True):
    x, values, diffs, colors = prepare_data(trim=trim)
    fig, ax = plt.subplots()
    ax.set_ylim(0, 1)
    for idx, (name, dy) in enumerate(diffs.items()):
        bottom = 0
        for i, segment in enumerate(dy):
            bar = ax.bar(name, segment, bottom=bottom, color=colors[i])
            bar_width = bar[0].get_width()
            # if i > 0:
            #     # Plot a dotted line connecting the corresponding segments
            #     ax.plot([idx - 1, idx], [bottom, bottom], linestyle=':', color='black')
            if idx > 0:
                # Get previous bar's corresponding segment
                prev_bottom = sum(list(diffs.values())[idx - 1][:i])
                # Plot a dotted line connecting the corresponding segments
                w = bar_width / 2
                x, y = [(idx-1)+w, idx-w], [prev_bottom, bottom]
                ax.plot(x,y, ls=':', c='black', lw=.8)
            bottom += segment
    # Set the tick positions and labels for the x-axis
    ax.set_xticks(range(len(diffs.keys())))
    # Set the rotation angle for the x-axis tick labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=60)


def plot_specific_function(name, trim=True):
    x, values, diffs, colors = prepare_data(name, trim)
    for name, y in values.items():
        plot(x, y, name)


# =============================================================================
# Auxilliary functions
# =============================================================================
def evaluate_function(x, func_name, trim=True, **kwargs):
    func = funcs.get(func_name)
    if func is None:
        raise ValueError(f"Function with name '{func_name}' not found")
    # Call the selected function
    trim_funcs = ('root', 'circle')
    if (func_name in trim_funcs) and trim:
        x = np.linspace(0, 1, len(x) + 1)
        y = func(x, **kwargs)[1:]
        y = (y-min(y))/(max(y)-min(y))
    else:
        y = func(x, **kwargs)
    return y


def prepare_data(specific_func=None, trim=True):
    # # Generate an array of equally spaced points between 0 and 1
    x = np.linspace(0, 1, 11)

    # Calculate function values for each function in the funcs dictionary
    if specific_func is None or specific_func == 'all':
        values = {name: evaluate_function(x, name, trim) for name in funcs.keys()}
    else:
        values = {specific_func: evaluate_function(x, specific_func, trim)}

    # Compute the differences between consecutive values for each function
    diffs = {k: np.diff(y) for k, y in values.items()}

    # # Get the first segment for each function
    # first_segments = {k: y[1] for k, y in values.items()}
    # # Sort the functions based on the first segment
    # sorted_list = sorted(diffs.items(), key=lambda item: first_segments[item[0]])
    # diffs = dict(sorted_list)
    # sorted_list = ['linear','poly','quadratic','cubic','easeOutExpo', 'circle','root']
    # diffs = {name: diffs[name] for name in sorted_list}

    # # First sort by y[1]
    # # sort1 = {k: y[1] for k, y in values.items()}
    # sort1 = {k: (y[1],y[2]) for k, y in values.items()}
    # sort1 = sorted(sort1.items(), key=lambda item: item[1])
    # # sort1 = sorted(values.items(), key=lambda item: sort1[item[0]])

    # for i,(k,v) in enumerate(sort1):
    #     v_next = sort1[i+1][1][1]
    #     v_curr = v[1]
    #     if v_curr > v_next:
    #         change_index = i+1
    #         break

    # sort2 = sort1[change_index:]
    # sort2 = {k: y[2] for k, y in sort2}
    # sort2 = sorted(sort2.items(), key=lambda item: item[0])
    # sorted_list = sort1[:change_index]
    # # sort2 = sorted(values.items(), key=lambda item: sort2[item[0]])

    # Create sort1 dictionary
    sort1 = {k: (y[1], y[2]) for k, y in values.items()}

    # Sort sort1 by y[1] in ascending order
    sort1_keys_sorted = sorted(sort1, key=lambda k: sort1[k][0])

    # Find the change_index
    for i, k in enumerate(sort1_keys_sorted[:-1]):
        curr_y2 = sort1[k][1]
        next_y2 = sort1[sort1_keys_sorted[i + 1]][1]
        if curr_y2 > next_y2:
            change_index = i + 1
            break

    # Divide the list into two parts
    first_part = sort1_keys_sorted[:change_index]
    second_part = sort1_keys_sorted[change_index:]

    # Sort the second part by y[2] in descending order
    second_part_sorted = sorted(second_part, key=lambda k: sort1[k][1], reverse=True)

    # Combine the two parts to get the final order
    final_order = first_part + second_part_sorted

    # Create a new dictionary sorted according to the final order of keys
    sorted_values = {k: values[k] for k in final_order}
    diffs = {k: diffs[k] for k in sorted_values.keys()}


    # # Find the index where the order changes from ascending to descending
    # change_index = None
    # for i in range(1, len(sorted_list)):
    #     if sort1[sorted_list[i][0]] < sort1[sorted_list[i - 1][0]]:
    #         change_index = i
    #         break

    # # Sort the remaining elements in descending order based on y[2]
    # if change_index is not None:
    #     sort2 = {k: y[2] for k, y in values.items()}
    #     sorted_list = (
    #         sorted_list[:change_index]
    #         + sorted(
    #             sorted_list[change_index:],
    #             key=lambda item: sort2[item[0]],
    #             reverse=True
    #         )
    #     )
    # diffs = dict(sorted_list)

    # Get the default color cycle
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # Generate additional colors using the viridis colormap
    additional_colors = plt.cm.viridis(np.linspace(0, 1, len(x)))
    additional_colors = [matplotlib.colors.to_hex(color) for color in additional_colors]
    # Combine default and additional colors
    colors = default_colors + additional_colors

    # Return the prepared data
    return x, values, diffs, colors


def plot(x,y,name):
    fig,ax = plt.subplots()
    ax.plot(x,y,'.-',label=name)
    for xx,yy in zip(x,y):
        ax.axhline(yy, xmin=0, xmax=xx, ls=':')
    ax.set_ylim(0,1)
    ax.set_xlim(0,1)
    ax.legend(loc='lower right')
    ax.set_aspect('equal', 'box')


# =============================================================================
# Call MAIN
# =============================================================================
if __name__=='__main__':
    import matplotlib.pyplot as plt
    import matplotlib.colors
    main()
