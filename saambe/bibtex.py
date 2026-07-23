# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Natl. Center of Biotechnology CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
@article{Pahari2020,
AUTHOR = {Pahari, Swagata and Li, Gen and Murthy, Adithya Krishna and Liang, Siqi and Fragoza, Robert and Yu, Haiyuan and Alexov, Emil},
TITLE = {SAAMBE-3D: Predicting Effect of Mutations on Protein–Protein Interactions},
JOURNAL = {International Journal of Molecular Sciences},
VOLUME = {21},
YEAR = {2020},
NUMBER = {7},
ARTICLE-NUMBER = {2563},
URL = {https://www.mdpi.com/1422-0067/21/7/2563},
PubMedID = {32272725},
ISSN = {1422-0067},
ABSTRACT = {Maintaining wild type protein–protein interactions is essential for the normal function of cell and any mutation that alter their characteristics can cause disease. Therefore, the ability to correctly and quickly predict the effect of amino acid mutations is crucial for understanding disease effects and to be able to carry out genome-wide studies. Here, we report a new development of the SAAMBE method, SAAMBE-3D, which is a machine learning-based approach, resulting in accurate predictions and is extremely fast. It achieves the Pearson correlation coefficient ranging from 0.78 to 0.82 depending on the training protocol in benchmarking five-fold validation test against the SKEMPI v2.0 database and outperforms currently existing algorithms on various blind-tests. Furthermore, optimized and tested via five-fold cross-validation on the Cornell University dataset, the SAAMBE-3D achieves AUC of 1.0 and 0.96 on a homo and hereto-dimer test datasets. Another important feature of SAAMBE-3D is that it is very fast, it takes less than a fraction of a second to complete a prediction. SAAMBE-3D is available as a web server and as well as a stand-alone code, the last one being another important feature allowing other researchers to directly download the code and run it on their local computer. Combined all together, SAAMBE-3D is an accurate and fast software applicable for genome-wide studies to assess the effect of amino acid mutations on protein–protein interactions. The webserver and the stand-alone codes (SAAMBE-3D for predicting the change of binding free energy and SAAMBE-3D-DN for predicting if the mutation is disruptive or non-disruptive) are available.},
DOI = {10.3390/ijms21072563}
}
"""