import logging
import platform
import test
from example import casestudy2
import figures_draw
import herb_herb_pairs
import permuation_herbs
from example import covid19

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

    # casestudy2.main()
    # figures_draw.main()
    # herb_herb_pairs.main()
    # permuation_herbs.main()
    # figures_draw.main()
    covid19.main()


