#include "../src/mesher.h"

#include <gtest/gtest.h>

#include <sstream>

#include "../src/vec2.h"

using namespace EllipticalPAN;

TEST(MesherTest, MesherTest1) {
    Mesher<Vec2<double>, double> Mesh(3, 3);
    Mesh.SetPoint(0, 0, Vec2<double>(0, 0));
    Mesh.SetPoint(1, 0, Vec2<double>(1, 0));
    Mesh.SetPoint(2, 0, Vec2<double>(1, 1));
    Mesh.SetPoint(2, 1, Vec2<double>(0.75, 1));
    Mesh.SetPoint(2, 2, Vec2<double>(0.5, 1));
    Mesh.SetPoint(1, 2, Vec2<double>(0.5, 0.5));
    Mesh.SetPoint(0, 2, Vec2<double>(0, 0.5));
    Mesh.SetPoint(0, 1, Vec2<double>(0, 0.25));
    Mesh.Generate();
    ASSERT_LE((Mesh.GetPoint(1, 1) - Vec2<double>(0.634615, 0.365385)).norm(),
              1e-5);
}

TEST(MesherTest, MesherTest2) {
    Mesher<Vec2<double>, double> Mesh(3, 3);
    Mesh.SetPoint(0, 0, Vec2<double>(0, 0));
    Mesh.SetPoint(1, 0, Vec2<double>(1, 0));
    Mesh.SetPoint(2, 0, Vec2<double>(2, 0));
    Mesh.SetPoint(2, 1, Vec2<double>(2, 2));
    Mesh.SetPoint(2, 2, Vec2<double>(2, 4));
    Mesh.SetPoint(1, 2, Vec2<double>(1, 3));
    Mesh.SetPoint(0, 2, Vec2<double>(0, 2));
    Mesh.SetPoint(0, 1, Vec2<double>(0, 1));
    Mesh.Generate();
    ASSERT_LE((Mesh.GetPoint(1, 1) - Vec2<double>(1, 1.39286)).norm(), 1e-5);
}
