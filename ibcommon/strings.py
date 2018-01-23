
import random as rnd
import string


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    """
    Create random string
    :param size: number of characters
    :param chars: characters to pick from
    :return:
    """
    return ''.join(rnd.choice(chars) for x in range(size))
