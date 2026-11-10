import plantbox as pb


def set_initial_heading(organ, heading):
    state = organ.getState()
    state["partialIHeading"] = list(heading)
    organ.setState(state)


def relative_nodes(organ, absolute_nodes):
    nodes = [pb.Vector3d(*absolute_nodes[0])]
    heading = organ.getInitialHeading()
    for previous, current in zip(absolute_nodes, absolute_nodes[1:]):
        delta = pb.Vector3d(
            current[0] - previous[0],
            current[1] - previous[1],
            current[2] - previous[2],
        )
        frame = pb.Matrix3d.ons(pb.Vector3d(heading))
        nodes.append(frame.inverse().times(delta))
        heading = delta
        heading.normalize()
    return nodes
