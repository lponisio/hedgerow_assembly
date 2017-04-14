
## if this were in R I would list all the files then divide by using
## grep on the names.... then loop over the name groups

import glob
from subprocess import call


pairs = glob.glob("*_.pairs")








## BACI
python runNetworkChangePoint.py graph-names.lut 4 'Barger_2006.pairs Barger_2007.pairs Barger_2008.pairs Barger_2009.pairs Barger_2011.pairs Barger_2012.pairs Barger_2013.pairs Barger_2014.pairs Barger_2015.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Butler_2006.pairs Butler_2007.pairs Butler_2008.pairs Butler_2009.pairs Butler_2011.pairs Butler_2012.pairs Butler_2013.pairs Butler_2014.pairs Butler_2015.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Hrdy_2007.pairs Hrdy_2008.pairs Hrdy_2009.pairs Hrdy_2011.pairs Hrdy_2012.pairs Hrdy_2013.pairs Hrdy_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'MullerB_2006.pairs MullerB_2007.pairs MullerB_2008.pairs MullerB_2009.pairs MullerB_2011.pairs MullerB_2012.pairs MullerB_2013.pairs MullerB_2014.pairs MullerB_2015.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Sperandio_2007.pairs Sperandio_2008.pairs Sperandio_2009.pairs Sperandio_2011.pairs Sperandio_2012.pairs Sperandio_2013.pairs Sperandio_2014.pairs' -p 'baci/'

## Controls
python runNetworkChangePoint.py graph-names.lut 4 'BC2_2007.pairs BC2_2008.pairs BC2_2009.pairs BC2_2011.pairs BC2_2012.pairs BC2_2013.pairs BC2_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Chamberlain_2007.pairs Chamberlain_2008.pairs Chamberlain_2011.pairs Chamberlain_2012.pairs Chamberlain_2013.pairs Chamberlain_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'DQU_2006.pairs DQU_2007.pairs DQU_2008.pairs DQU_2009.pairs DQU_2011.pairs DQU_2012.pairs DQU_2013.pairs DQU_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Turkovich_2006.pairs Turkovich_2007.pairs Turkovich_2008.pairs Turkovich_2009.pairs Turkovich_2011.pairs Turkovich_2012.pairs Turkovich_2013.pairs Turkovich_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'USS_2006.pairs USS_2007.pairs USS_2008.pairs USS_2009.pairs USS_2011.pairs USS_2012.pairs USS_2013.pairs USS_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'H16_2006.pairs H16_2007.pairs H16_2008.pairs H16_2009.pairs H16_2011.pairs H16_2013.pairs H16_2014.pairs' -p 'baci/'

##
python runNetworkChangePoint.py graph-names.lut 4 'Gregory_2007.pairs Gregory_2009.pairs Gregory_2011.pairs Gregory_2012.pairs Gregory_2013.pairs Gregory_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Spiva_2007.pairs Spiva_2008.pairs Spiva_2011.pairs Spiva_2012.pairs Spiva_2013.pairs Spiva_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'MC1_2007.pairs MC1_2008.pairs
MC1_2009.pairs MC1_2011.pairs MC1_2012.pairs MC1_2013.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Joe_2009.pairs Joe_2010.pairs Joe_2012.pairs Joe_2013.pairs Joe_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Tractor_2009.pairs Tractor_2011.pairs Tractor_2012.pairs Tractor_2013.pairs Tractor_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Rominger_2009.pairs Rominger_2010.pairs Rominger_2011.pairs Rominger_2012.pairs Rominger_2013.pairs Rominger_2014.pairs' -p 'baci/'

## mature
python runNetworkChangePoint.py graph-names.lut 4 'Fong_2008.pairs Fong_2009.pairs Fong_2010.pairs Fong_2011.pairs Fong_2012.pairs Fong_2013.pairs Fong_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'Harlan_2009.pairs Harlan_2011.pairs Harlan_2012.pairs Harlan_2013.pairs Harlan_2014.pairs' -p 'baci/'

python runNetworkChangePoint.py graph-names.lut 4 'MullerM_2010.pairs MullerM_2011.pairs MullerM_2012.pairs MullerM_2013.pairs MullerM_2014.pairs' -p 'baci/'

