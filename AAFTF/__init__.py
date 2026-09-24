"""Init placeholder."""

import logging

# Library code only logs; AAFTF_main.main() configures handlers via setup_logging().
logging.getLogger(__name__).addHandler(logging.NullHandler())
