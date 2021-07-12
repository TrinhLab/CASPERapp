import logging

"""File stores all global variables, mainly used for file information/directory monitoring"""

# setup logger in global space
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

OPERATING_SYSTEM_ID = ""

CSPR_DB = ""

appdir = ""