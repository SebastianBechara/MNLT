#!/usr/bin/env python3

#Author: Estienne C. Swart
#Run after_ParTIES.py with the --help switch for usage information

from sys import argv
from matplotlib import *
from numpy import mean, std, array, arange, linspace
import matplotlib.cm as cm
import matplotlib.style as style
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LogNorm
from scipy import stats
from scipy.odr import ODR, Model, RealData
from statsmodels.nonparametric.smoothers_lowess import lowess
import argparse, sys

def argParse():
  parser = argparse.ArgumentParser(description=
"""Using the IES retention scores produced by ParTIES, generate a correlation matrix.
As an example, the following command should generate an image like the one in
Fig 3A of Swart et al. (all text between quotes type on a single line):
'python3.6 after_ParTIES.py --input_file ies_retention.tab
--output_matrix_image ies_retention_knockdown_matrix.png
--experiment_columns 1 2 3 10 6 5 7 4 9 12 --control_columns 11 12 13 14'
""", add_help=False,usage=argparse.SUPPRESS,formatter_class=argparse.RawTextHelpFormatter)

  required_arg_group = parser.add_argument_group("required options")
  required_arg_group.add_argument("--input_file", 
    help="Tab delimited file of IES retention scores from ParTIES.") #required=True,
  required_arg_group.add_argument("--output_matrix_image", 
    help="Output file for correlation matrix image. The suffix of the filename \
will determine the image type (e.g. .png = png file). For a few experiments pdf \
images can be generated, but for many onscreen rendering becomes too slow, \
and so either jpg or png files should be created.") #required=True,
  required_arg_group.add_argument('--experiment_columns', metavar='E', type=int, nargs='+',
     help='One or more columns of IES retention scores from an \
experimental series of knockdowns. [zero based index]') #required=True,
  required_arg_group.add_argument('--control_columns', metavar='N', type=int, nargs='+',
     help='One or more columns of IES retention scores from \
negative controls of experimental knockdowns. [zero based index]')#required=True,

  optional_arg_group = parser.add_argument_group("Optional options")
  optional_arg_group.add_argument("--use_pearson", action="store_true",
    help="Add this option if Pearson's correlation is to be used instead of \
Spearman's correlation. In this case no p-value testing is performed.")
  optional_arg_group.add_argument("--alpha", default=0.01, type=float,
    help="Alpha used for p-value testing of Spearman's correlation coefficient \
(default: %(default)s); a Bonferroni correction is applied to this before \
testing, and when p > alpha/(number of hypotheses), a '^' character is shown \
next to the calculated correlation coefficient in the upper diagonal matrix. \
For example, for 10 experimental knockdowns the number of hypotheses is 45.")
  optional_arg_group.add_argument("--deactivate_ODR", action="store_true", default=False,
    help="Use this flag to turn of ODR for a speedup (default: %(default)s).")
  optional_arg_group.add_argument("--image_resolution", default=400, type=int,
    help="Change image resolution (default (in DPI): %(default)s).")
  optional_arg_group.add_argument("--base_font_size", default=14, type=int,
    help="Change font size (default: %(default)s).")
  optional_arg_group.add_argument("--font_style", default='italic',
    help="Change font style (default: %(default)s); options = normal, italic \
or oblique.")
  optional_arg_group.add_argument("--regression_line_width", default=1.2, type=int,
    help="Thickness of regression lines (default: %(default)s).")
  optional_arg_group.add_argument("--show_colorbar",
    help="Show colorbar scale of hexagonal bins.", action="store_true")
  optional_arg_group.add_argument("--show_axis_labels",
    help="Show labels on axes", default=True)
  optional_arg_group.add_argument("--histogram_max", default=9000, type=int,
    help="Maximum value of y-axis of the matrix diagonal's histograms (default: %(default)s).")
  optional_arg_group.add_argument("--hexbin_vmin", default=1, type=int,
    help="Minimum value of hexagonal bin counts (log10 scale) (default: %(default)s).")
  optional_arg_group.add_argument("--hexbin_vmax", default=1000, type=int,
    help="Maximum value of hexagonal bin counts (log10 scale) (default: %(default)s).")

  logging_arg_group = parser.add_argument_group("MISC options")
  logging_arg_group.add_argument("--verbose", help="Display ODR output \
that may be useful for diagnostic purposes.", action="store_true")
  logging_arg_group.add_argument("--help", help="Display this message \
and exit.", action="help")



  return parser

def main(args=None):
  parser= argParse()
  num=1
  if args:
    num=2
  if not args:
    args = parser.parse_args()
  if len(sys.argv) == num:
    parser.print_help(sys.stderr)
    sys.exit()

  style.use('classic')
  #for the moment revert to classic style due to problems with axis
  #ticking in matplotlib 2.0.0
  rcParams['font.sans-serif'] = 'Arial'
  tick_font_size = args.base_font_size - len(args.experiment_columns)
  label_font_size = args.base_font_size - 0.75*len(args.experiment_columns)

  lines = open(args.input_file).readlines()
  data_lines = lines[1:]
  experiment_names = lines[0].split('\t')
  experiment_index = args.experiment_columns
  experiments = [experiment_names[ind] for ind in experiment_index]
  n_subfigs = len(experiments)

  corrected_p_threshold = args.alpha/sum(range(len(experiments)))
  #Bonferroni correction

  def fit_func(B, x):
    """y = m*x + b"""
    return B[0]*x + B[1]

  linear = Model(fit_func)

  fig = plt.figure()
  color_generator = iter(cm.magma(arange(0, 1, 1.0/n_subfigs)))

  z = 0
  for i in range(n_subfigs):
    experiment1 = experiments[i]
    for j in range(n_subfigs):
      experiment2 = experiments[j]
      z += 1

      vals_1 = {}
      for line in data_lines:
        atoms = line.split()
        if atoms[experiment_index[i]] != 'NA':
          #exclude NA values generated by ParTIES
          score = float(atoms[experiment_index[i]])
          ctrl_mean = mean([float(atoms[n]) for n in range(len(atoms)) if n in args.control_columns])
          vals_1[atoms[0]] = max(score - ctrl_mean, 0)

      vals_2 = {}
      for line in data_lines:
        atoms = line.split()
        if atoms[experiment_index[j]] != 'NA':
          #exclude NA values generated by ParTIES
          score = float(atoms[experiment_index[j]])
          ctrl_mean = mean([float(atoms[n]) for n in range(len(atoms)) if n in args.control_columns])
          vals_2[atoms[0]] = max(score - ctrl_mean, 0)

      common_keys = set(vals_1.keys()) & set(vals_2.keys())

      x, y = [], []
      for akey in common_keys:
        x.append(vals_2[akey])
        y.append(vals_1[akey])

      slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
      rho, rho_p = stats.spearmanr(x, y)

      if i > j:
        if len(x) > 0 and len(y) > 0:
          plt.subplot(n_subfigs, n_subfigs, z)
          plt.hexbin(x, y, norm=LogNorm(vmin=args.hexbin_vmin, vmax=args.hexbin_vmax), mincnt=1, cmap=cm.viridis)

          max_max = max(max(x), max(y))
          xi = arange(0, max_max, 0.01)
          line = slope*xi+intercept

          mydata = RealData(array(x), array(y), sx=std(x), sy=std(y))

          if not args.deactivate_ODR:
            myodr = ODR(mydata, linear, beta0=[slope, intercept])
            ODR_out = myodr.run()

            if args.verbose:
              ODR_out.pprint()

            ODR_fitted_line = fit_func(ODR_out.beta, array(x))

          lowess_arr = lowess(y, x, delta=0.01, return_sorted=True)
          ys = lowess_arr[:,1]
          xs = lowess_arr[:,0]

          plt.plot(xs, ys,'orange',lw=args.regression_line_width)
          plt.plot(xi, line, ls='-', color='red', lw=args.regression_line_width)
          if not args.deactivate_ODR:
            plt.plot(x, ODR_fitted_line, "--", color='gray', lw=args.regression_line_width)

          plt.axis(xmax=1.01, ymax=1.01, ymin=-0.01, xmin=-0.01)

          ax = plt.gca()
          current_xticklabels = ax.get_xticks()
          fixed_xticklabels = ["%.1f" % (round(current_label, 1)) for current_label in current_xticklabels]
          current_yticklabels = ax.get_yticks()
          fixed_yticklabels = ["%.1f" % (round(current_label, 1)) for current_label in current_yticklabels]

          ax.yaxis.set_major_locator(mticker.FixedLocator(list(current_yticklabels)))
          ax.xaxis.set_major_locator(mticker.FixedLocator(list(current_xticklabels)))

          ax.set_xticklabels([])
          ax.set_yticklabels([])


          plt.tick_params(
          axis='both',
          which='both',
          top=False,
          right=False)

          if j == 0:
            #ax.set_xticklabels(fixed_xticklabels, fontsize=tick_font_size)
            ax.set_yticklabels(fixed_yticklabels, fontsize=tick_font_size)

            if args.show_axis_labels:
              ax.set_ylabel(experiment1, fontsize=label_font_size,
                rotation=90, style=args.font_style)
          if i == n_subfigs -1:
            ax.set_xticklabels(fixed_xticklabels, fontsize=tick_font_size, rotation=45)
            if args.show_axis_labels:
              ax.set_xlabel(experiment2, fontsize=label_font_size, style=args.font_style)


      elif i == j:
        plt.subplot(n_subfigs, n_subfigs, z)

        plt.tick_params(
        axis='both',
        which='both',
        top=False,
        right=False)

        if len(x) > 0:
          plt.hist(x, bins=arange(0, 1, 0.025), linewidth=0.3, ec="white",
            fc=next(color_generator))
          plt.axis(ymin=0, ymax=args.histogram_max)

          ax = plt.gca()
          ax.set_yticks(ax.get_yticks()[::2])

          ax.set_xticks(arange(0, 1.2, 0.2))
          current_xticklabels = ax.get_xticks()
          fixed_xticklabels = ["%.1f" % (round(current_label, 1)) for current_label in current_xticklabels]
          current_yticklabels = ax.get_yticks()
          fixed_yticklabels = ["%.1f" % (round(current_label, 1)) for current_label in current_yticklabels]

          ax.yaxis.set_major_locator(mticker.FixedLocator(list(current_yticklabels)))
          ax.xaxis.set_major_locator(mticker.FixedLocator(list(current_xticklabels)))

          ax.set_xticklabels([], fontsize=tick_font_size)
          ax.set_yticklabels(fixed_yticklabels, fontsize=tick_font_size)

          if i == n_subfigs -1:
            ax.set_xticklabels(fixed_xticklabels, fontsize=tick_font_size, rotation=45)

          if i > 0:
            plt.tick_params(
            axis='both',
            labelleft=False)

      else:
        plt.subplot(n_subfigs, n_subfigs, z)
        ax = plt.gca()
        current_xticklabels = ax.get_xticks()
        fixed_xticklabels = ["%.1f" % (round(current_label, 1)) for current_label in current_xticklabels]

        ax.yaxis.set_major_locator(mticker.FixedLocator(list(current_yticklabels)))
        ax.xaxis.set_major_locator(mticker.FixedLocator(list(current_xticklabels)))
        ax.set_xticklabels(fixed_xticklabels, fontsize=tick_font_size)


        if not args.use_pearson:
          ax.set_facecolor(cm.YlGnBu(rho))
          if rho > 0.5:
            if rho_p < corrected_p_threshold:
              plt.text(0.2, 0.4, "%.2f" % round(rho, 3), fontsize=1.5*label_font_size, color='white')
            else:
              plt.text(0.2, 0.4, "%.2f^" % round(rho, 3), fontsize=1.5*label_font_size, color='white')
          else:
            if rho_p < corrected_p_threshold:
              plt.text(0.2, 0.4, "%.2f" % round(rho, 3), fontsize=1.5*label_font_size, color='black')
            else:
              plt.text(0.2, 0.4, "%.2f^" % round(rho, 3), fontsize=1.5*label_font_size, color='black')
        else:
          ax.set_facecolor(cm.YlGnBu(r_value))
          if r_value > 0.5:
            plt.text(0.2, 0.4, "%.2f" % round(r_value, 3), fontsize=1.5*label_font_size, color='white')
          else:
            plt.text(0.2, 0.4, "%.2f" % round(r_value, 3), fontsize=1.5*label_font_size, color='black')

        plt.tick_params(
        axis='both',
        which='both',
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labelleft=False,
        labelright=False)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

  if args.show_colorbar:
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.4, 0.015, 0.4])

    cb_label = 'correlation plot hexbin scale (log10)'
    cb = plt.colorbar(cax=cbar_ax, drawedges=False,
      ticks=range(args.hexbin_vmin, args.hexbin_vmax+1), label=cb_label)
    #current_yticklabels = cbar_ax.get_yticks()
    #fixed_yticklabels = ["%.1f" % (round(current_label, 1)) for current_label in current_yticklabels]

    cbar_ax.set_yticklabels(range(args.hexbin_vmin, args.hexbin_vmax+1), fontsize=tick_font_size)
    cbar_ax.set_ylabel(cb_label, fontsize=label_font_size)
    cb.outline.set_linewidth(1)
    #cb.ax.get_children()[4].set_linewidths(1)

  plt.savefig(args.output_matrix_image, dpi=((args.image_resolution)))


if __name__ == "__main__":
  main()



