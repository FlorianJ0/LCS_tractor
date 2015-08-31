__author__ = 'p0054421'


# see http://danishmujeeb.com/blog/2014/01/basic-sorting-algorithms-implemented-in-python
def insertion_sort(items):
    """ Implementation of insertion sort 
    :type items: list
    """
    for i in range(1, len(items)):
        j = i
        while j > 0 and items[j] < items[j - 1]:
            items[j], items[j - 1] = items[j - 1], items[j]
            j -= 1

    return items
