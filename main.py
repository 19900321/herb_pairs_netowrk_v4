import logging
import platform
import test
from example import casestudy2
import figures_draw
import herb_herb_pairs
import permuation_herbs

logger = logging.getLogger()
if platform.node() == 'calculon':
    logger.setLevel('INFO')
else:
    logger.setLevel('DEBUG')
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler()
ch.setLevel('INFO')
ch.setFormatter(formatter)
logger.addHandler(ch)


if __name__ == '__main__':
    # test.main()

    casestudy2.main()
    # figures_draw.main()
    # herb_herb_pairs.main()
    # permuation_herbs.main()
    # figures_draw.main()


# test
import pandas as pd
import numpy as np
pd.read_excel()