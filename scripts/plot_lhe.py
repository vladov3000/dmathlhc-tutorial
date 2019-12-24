import math
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt


class LHEParticle(object):
    """
    Class copied from https://github.com/lukasheinrich/pylhe/blob/master/src/pylhe/__init__.py
    """

    fieldnames = [
        "id",
        "status",
        "mother1",
        "mother2",
        "color1",
        "color2",
        "px",
        "py",
        "pz",
        "e",
        "m",
        "lifetime",
        "spin",
    ]

    def __init__(self, **kwargs):
        if not set(kwargs.keys()) == set(self.fieldnames):
            raise RuntimeError
        for k, v in kwargs.items():
            setattr(self, k, v)

    @classmethod
    def fromstring(cls, string):
        obj = cls(**dict(zip(cls.fieldnames, map(float, string.split()))))
        return obj

def transverse_momentum(particles):
    return math.sqrt(sum([sum([getattr(p, mu) for p in particles]) ** 2 for mu in ['px', 'py']]))


def invariant_mass(particles):
    E = sum([p.e for p in particles])
    p = math.sqrt(sum([sum([getattr(p, mu) for p in particles]) ** 2 for mu in ['px', 'py', 'pz']]))
    return math.sqrt(E ** 2 - p ** 2)


def find_dm(particles):
    out = []
    for p in particles:
        if 51 <= abs(p.id // (10 ** 5)) <= 53:
            out.append(p)
    return out


def get_event_particles(lhe_filepath):
    LHEevents = []
    for event, element in ET.iterparse(lhe_filepath):
        if element.tag == "event":
            data = element.text.split("\n")[1:-1]
            event_data, particles = data[0], data[1:]
            particle_objs = []
            for p in particles:
                particle_objs.append(LHEParticle.fromstring(p))
            LHEevents.append(particle_objs)
    return LHEevents


def plot_func(events, func):
    values = list(map(func, map(find_dm, events)))
    plt.title(func.__name__ + " Distribution")
    plt.xlabel(func.__name__)
    plt.ylabel('Events')
    plt.hist(values)
    plt.show()


def test():
    test_file = "../test_data/unweighted_events.lhe"
    events = get_event_particles(test_file)

    plot_func(events, transverse_momentum)
    plot_func(events, invariant_mass)


if __name__ == "__main__":
    test()
