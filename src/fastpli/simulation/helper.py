from .generation import LayerProperty


def TupleList2Layer(tuple_lists):
    fbp = []
    for tuple_list in tuple_lists:
        fbp.append([])
        for elm in tuple_list:
            fbp[-1].append(LayerProperty(elm[0], elm[1], elm[2], elm[3]))
    return fbp
