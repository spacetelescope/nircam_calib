python nircam_full_system_throughput_step_by_step_plots.py filteronly_modA.list --opticsfile NIRCam_optics_transmission_29Oct2015.csv -m a
python nircam_full_system_throughput_step_by_step_plots.py filteronly_modB.list --opticsfile NIRCam_optics_transmission_29Oct2015.csv -m b

python mod_ab_means.py

python multi_filter_throughput_plot.py nrc_and_ote_meanmod.list --nrc_optics --ote --module mean --ind_plots
python multi_filter_throughput_plot.py nrc_only_meanmod.list --nrc_optics --module mean --ind_plots
python multi_filter_throughput_plot.py filter_only_meanmod.list --module mean --ind_plots --filteronly
