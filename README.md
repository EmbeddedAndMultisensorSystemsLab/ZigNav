This is a DEMO of ZigNav, a wireless positioning system developed in the Embedded Multisensor Systems (EMS) research laboratory in Carleton University. ZigNav dynamically generate radiomaps to localize people/objects using received signal strength measurements. ZigNav uses cooperative measurements exchanged between wireless nodes to construct a dense radio map and to use the generated maps to localize mobile objects. ZigNav can work on any existing wireless infrastructure with minimal overhead. EMS has developed a hybrid approach that combines gaussian process regression models and pathloss models to generate online radiomaps from the cooperative measurements exchanged between wireless nodes. 

In this DEMO, the positioing in an indoor area with few wireless ZigBee nodes is demonstrated. Technically, the input logs to the system is only the cooperative measurements exchanged between wireless nodes. ZigNav uses this exchanged messages to dynamically generate radiomaps using a hybrid models approach. However, in this DEMO, the radiomaps have already been created. Therefore, the radiomap generation code is not given in this DEMO package. The DEMO show the localization process using k-nearest-neighbourhood (knn) algorithm. This is a 20-minutes stationary testing scenario. Dynamic scenarios and integration with other sensors such as inertial measurement units and vision sensors are currently under development.

Three folders are given

1) dynamic_radiomaps
	- hybrid_radiomaps --this folder contains the generated radiomaps using EMS hybrid approach
	- pathloss_model_radiomaps --this folder contains the generated radiomaps using traditional pathloss models approach
	- zero_mean_gpr__radiomaps --this folder contains the generated radiomaps using gaussian-process regression approach

2) experiment_data -- this folder contains information files about the testing area such as dimension, floor map, ground-truth wireless nodes pixel locations, and the recorded cooperative measurements exchanged between wireless nodes.

3) wireless_positioning --this folder contains the positioning code

For the hybrid radiomap generation method, check the following paper
M.M. Atia, M. Korenberg, A. Noureldin, “Dynamic Online-Calibrated Radio Maps for Indoor Positioning In Wireless Local Area Networks”, IEEE Transactions on Mobile Computing, Vol.12, no.9, pp. 1774 - 1787, 2013

%   Copyright (C) 2018, Mohamed Atia, all rights reserved.
%   The software is given under GNU Lesser General Public License
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this program. If not, see
%   <http://www.gnu.org/licenses/>.

%   For commercial use, additional features, or support, please contact mohamed.atia@carleton.ca
