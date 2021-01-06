from enum import Enum


# TODO LOGIC: Combine modes and action states into one class
class ActionState:
    '''
        ActionState defines the different states the possible actions (filter, optimize etc.) can be in
    '''
    # Using the Pymol convention for numbers
    # If I use 4 for CURRENTLY_DISABLED and 1 for ALWAZS DISABLED, PYMOL WILL ONLY USE THE GREEN FONT FOR BOTH
    CURRENTLY_ENABLED = 3
    CURRENTLY_DISABLED = 1
    ALWAYS_ENABLED = 2
    ALWAYS_DISABLED = 4
    # Only usefull during debugging


class EventType(Enum):
    '''
        EventType defines the different events that a selection can deal with
    '''

    SELECT = 'select'
    MOVE = 'move'
    SIZE = 'size'
    CONFIRM = 'confirm'
