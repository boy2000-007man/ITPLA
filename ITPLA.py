#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import division
import matplotlib.path as pt
import Box2D.b2 as b2
import collections as co
import operator as op
import random as rd
import numpy as np
import heapq as hq


def rand(polygon):
    border = pt.Path(polygon)
    x, y = zip(*polygon)
    while True:
        p = (rd.uniform(min(x), max(x)), rd.uniform(min(y), max(y)))
        if border.contains_point(p):
            return p


# a polygon area calculator (CCW)
# Shoelace formula
def area(x, y): return (np.roll(x, 1).dot(y) - np.roll(y, 1).dot(x)) / 2


def distance2line(p, l):
    return 2 * area(*zip(p, *l)) / np.linalg.norm(np.subtract(*l))


def body2worldvs(b): return map(
    lambda v: b.transform * v, b.fixtures[0].shape.vertices)


def vs2xy(vs, x, y): return sorted(map(lambda v: [distance2line(
    v, y), distance2line(v, x)], vs), key=op.itemgetter(0))


def x2y(x, xy): return np.array(x).dot(np.linalg.lstsq(
    *zip(*map(lambda (x, y): [[1, x], y], xy)))[0])


def surface_distance(vs, es, x, y):
    surf_dist = []
    for v in vs:
        for e in es:
            xy = vs2xy(e, (v, v + x), (v, v + y))
            if xy[0][0] < 0 < xy[1][0]:
                surf_dist.append(x2y([1, 0], xy))
    return surf_dist


def weight(surf_dist, dist): return surf_dist ** 4 + dist ** -4


def coeffient(surf_dist, dist): return (1 - surf_dist / dist) ** 2 - 1


def included_angle(*ls):
    return (np.arctan2(*ls[0]) - np.arctan2(*ls[1]) + np.pi / 2) % np.pi - np.pi / 2

CCW90 = np.array([[0, -1], [1, 0]])

NEIGHBOUR_FLAG = True
OVERLAP_FLAG = None

lam = 0

TARGET_FPS = 24
TIME_STEP = 1 / TARGET_FPS


class ITPLA:
    def __init__(self, polygon, module, B_P):
        self.world = b2.world(gravity=(0, 0), doSleep=False)
        self.borders = map(lambda v: self.world.CreateStaticBody(shapes=b2.edgeShape(vertices=(
            v - np.mean(v, axis=0)).tolist()), position=np.mean(v, axis=0).tolist()), zip(np.roll(polygon, 2), polygon))
        self.modules = map(lambda _: self.world.CreateDynamicBody(shapes=b2.polygonShape(
            vertices=module.tolist()), position=rand(polygon), angle=rd.random() * np.pi), range(B_P))
        for b in self.borders + self.modules:
            b.fixtures[0].sensor = True

    def __iter__(self):
        return self

    def next(self):
        self.step()
        return self.world

    def step(self):
        num = len(self.modules)
        position = map(lambda m: m.position, self.modules)
        vertices = map(body2worldvs, self.modules)
        distance = map(lambda vs: map(lambda p: map(
            lambda v: np.linalg.norm(p - v), vs), position), vertices)
        edge = map(lambda ds: map(lambda d: tuple(sorted(zip(
            *hq.nsmallest(2, zip(range(len(d)), d), key=op.itemgetter(-1)))[0])), ds), distance)
        neighbour = map(lambda _: co.defaultdict(list), range(num))
        for i in range(num):
            for j in range(num):
                if i != j:
                    neighbour[i][edge[i][j]].append(
                        (j, np.linalg.norm(position[i] - position[j])))
        neighbour = map(lambda i: dict(map(lambda e: [e, min(
            neighbour[i][e], key=op.itemgetter(-1))[0]], neighbour[i])), range(num))

        for i in range(num):
            mi, pi, ivs = self.modules[i], position[i], vertices[i]
            ies = zip(np.roll(ivs, 2), ivs)

            weights, linearVelocity, angularVelocity = [], [], []

            for ei, j in neighbour[i].iteritems():
                mj, pj, jvs = self.modules[j], position[j], vertices[j]
                jes = zip(np.roll(jvs, 2), jvs)

                dist = np.linalg.norm(pi - pj)
                r = np.array(pi - pj)
                r /= np.linalg.norm(r)
                n = CCW90.dot(r)
                surf_dist = min(surface_distance(ivs, jes, n, r) +
                                surface_distance(jvs, ies, -n, r))
                k_r = coeffient(surf_dist, dist)

                ej = edge[j][i]
                eivs = op.itemgetter(*ei)(vertices[i])
                ejvs = op.itemgetter(*ej)(vertices[j])
                k_n = distance2line(np.mean(eivs, axis=0), ((
                    0, 0), r)) - distance2line(np.mean(ejvs, axis=0), ((0, 0), r))
                if NEIGHBOUR_FLAG:
                    weights.append(weight(surf_dist, dist))
                    linearVelocity.append(k_r * r)
                    linearVelocity[-1] -= k_n * n
                    ia = included_angle(np.subtract(*eivs), np.subtract(*ejvs))
                    angularVelocity.append(ia)

            for c in filter(lambda c: c.contact.touching, mi.contacts):
                ovs = body2worldvs(c.other)
                if isinstance(c.other.fixtures[0].shape, b2.edgeShape):
                    dist = distance2line(pi, ovs)
                    v = np.subtract(*ovs)
                    v /= np.linalg.norm(v)
                    n = CCW90.dot(v)
                    if 0 < dist:
                        surf_dist = min(
                            [0] + surface_distance(ivs, [ovs], v, n) + surface_distance(ovs, ies, -v, n))
                        k_n = coeffient(surf_dist, dist)
                        if OVERLAP_FLAG:
                            weights.append(weight(surf_dist, dist))
                            linearVelocity.append(-k_n * n)
                            angularVelocity.append(included_angle(np.subtract(
                                *ies[np.argmin(map(lambda e: distance2line(np.mean(e, axis=0), ovs), ies))]), np.subtract(*ovs)))
                elif isinstance(c.other.fixtures[0].shape, b2.polygonShape):
                    pass
            mi.linearVelocity *= lam
            mi.angularVelocity *= lam
            if weights:
                weights = np.array(weights) / sum(weights)
                mi.linearVelocity += np.average(linearVelocity,
                                                axis=0, weights=weights) / 10
                mi.angularVelocity += np.average(angularVelocity,
                                                 axis=0, weights=weights) / 10

        self.world.Step(timeStep=TIME_STEP, velocityIterations=6,
            positionIterations=2)

__all__ = ('ITPLA', 'area')