# MMV-AMP-with-Side-Information
This code is for paper: [On Massive IoT Connectivity with Temporally-Correlated User Activity](https://arxiv.org/pdf/2101.11344.pdf).
# Abstract
This paper considers joint device activity detection and channel estimation in Internet of Things (IoT) networks, where a large number of IoT devices exist but merely a random subset of them become active for short-packet transmission at each time slot. In particular, we propose to leverage the temporal correlation in user activity, i.e., a device active at the previous time slot is more likely to be still active at the current moment, to improve the detection performance. Despite the temporally-correlated user activity in consecutive time slots, it is challenging to unveil the connection between the activity pattern estimated previously, which is imperfect but the only available side information (SI), and the true activity pattern at the current moment due to the unknown estimation error. In this work, we manage to tackle this challenge under the framework of approximate message passing (AMP). Specifically, thanks to the state evolution, the correlation between the activity pattern estimated by AMP at the previous time slot and the real activity pattern at the previous and current moment is quantified explicitly. Based on the well-defined temporal correlation, we further manage to embed this useful SI into the design of the minimum mean-squared error (MMSE) denoisers and loglikelihood ratio (LLR) test based activity detectors under the AMP framework. Theoretical comparison between the SI-aided AMP algorithm and its counterpart without utilizing temporal correlation is provided. Moreover, numerical results are given to show the significant gain in activity detection accuracy brought by the SI-aided algorithm.
# Read Me
model.m is the main function <br>
noisyCAMPmmseforKLS.m and noisyCAMPmmseforKLSwithsi.m are used to calculate \eta for the case without si and with si respectively <br>
To plot the tradeoff curve, you need to tune on the threshold case by case. If you have any question or find any bug, you can contact me at: qipeng.wang@connect.polyu.hk 
# Citation
@article{wang2021massive,<br>
  title={On Massive IoT Connectivity with Temporally-Correlated User Activity},<br>
  author={Wang, Qipeng and Liu, Liang and Zhang, Shuowen and Lau, Francis},<br>
  journal={arXiv preprint arXiv:2101.11344},<br>
  year={2021}<br>
}
