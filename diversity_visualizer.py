import os
import sys
import shutil
import logging
import imp

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
py_src = os.path.join(home_directory, "py/pipeline/")
visualizer_dir = os.path.join(home_directory, "py/diversity_stats_visualizer")

sys.path.append(py_src)
sys.path.append(visualizer_dir)

import process_cfg
import support
import utils
import visualize_vj_stats
import visualize_cdr_stats
import visualize_shm_stats
import html_report_writer

###############################################################################
class VisualizerConfig:
    # input info
    tool_name = 'diversity_visualizer'
    cdr_details = "cdr_details.txt"
    shm_details = "shm_details.txt"
    # output files
    visualizer_dir = "visualization"
    plot_dir = "plots"
    html_report = "annotation_report.html"
    # packages
    matplotlib_package = 'matplotlib'
    seaborn_package = 'seaborn'
    pandas_package = 'pandas'

def CheckInputDirFatal(input_dir):
    if not os.path.isdir(input_dir):
        print "ERROR: " + input_dir + " is not a directory"
        sys.exit(1)

def CreateOutputDir(divan_dir):
    vis_dir = os.path.join(divan_dir, VisualizerConfig.visualizer_dir)
    if not os.path.exists(vis_dir):
        os.mkdir(vis_dir)
    return vis_dir

def CreateLog(output_dir):
    log = logging.getLogger(VisualizerConfig.tool_name)
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    # adding log to file
    log_filename = os.path.join(output_dir, VisualizerConfig.tool_name + '.log')
    if os.path.exists(log_filename):
        os.remove(log_filename)
    log_handler = logging.FileHandler(log_filename, mode='a')
    log.addHandler(log_handler)
    return log

def CheckInputFilesFatal(input_dir, log):
    cdr_details = os.path.join(input_dir, VisualizerConfig.cdr_details)
    if not os.path.exists(cdr_details):
        log.info("ERROR: CDR annotation file " + cdr_details + " was not found")
        sys.exit(1)
    shm_details = os.path.join(input_dir, VisualizerConfig.shm_details)
    if not os.path.exists(shm_details):
        log.info("ERROR: SHM annotation file " + shm_details + " was not found")
        sys.exit(1)
    return cdr_details, shm_details

def CheckModuleExistanceFatal(log):
    utils.CheckPackageFatal(VisualizerConfig.matplotlib_package, log)
    utils.CheckPackageFatal(VisualizerConfig.pandas_package, log)
    utils.CheckPackageFatal(VisualizerConfig.seaborn_package, log)

#######################################################################################
def main(input_dir, output_dir, output_log):
    CheckModuleExistanceFatal(output_log)
    cdr_details, shm_details = CheckInputFilesFatal(input_dir, output_log)
    plot_dir = os.path.join(output_dir, VisualizerConfig.plot_dir)
    os.mkdir(plot_dir)
    output_log.info("==== Visualization of diversity statistics ====")
    output_log.info("\n== Output VJ statistics ==")
    visualize_vj_stats.main(["", cdr_details, shm_details, output_dir, output_log])

    output_log.info("\n== Output CDR / FR statistics ==")
    visualize_cdr_stats.main(cdr_details, plot_dir, output_log)

    output_log.info("\n== Output SHM statistics ==")
    visualize_shm_stats.main(shm_details, plot_dir, output_dir, output_log)

    output_log.info("\n==== Annotation report creation ====")
    html_fname = os.path.join(output_dir, VisualizerConfig.html_report)
    html_report_writer.main(cdr_details, shm_details, plot_dir, html_fname, output_log)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Invalid input arguments"
        print "python diversity_visualizer.py diversity_analyzer_dir"
        sys.exit(1)
    diversity_analyzer_dir = sys.argv[1]
    CheckInputDirFatal(diversity_analyzer_dir)
    output_dir = CreateOutputDir(diversity_analyzer_dir)
    output_log = CreateLog(output_dir)
    main(diversity_analyzer_dir, output_dir, output_log)