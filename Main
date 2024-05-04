#include <embree3/rtcore.h>
#include <embree3/rtcore_ray.h>

#include "integra.h"

#include "base/str.hpp"
#include "base/file.hpp"
#include "base/marray.hpp"
#include "math/matrix43.hpp"
#include "math/vect2.hpp"
#include "math/vect3.hpp"
#include "base/matrix.hpp"

START_C_DECLS
#include "ievl.h"
// #include "iosl.h"
#include "icol.h"
#include "imal.h"
#include "iifl.h"
#include "itoliifl.h"
#include "itoliiff.h"
#include <math/rnd.hpp>
// #include "itolscan.h"
// #include "suffix.h"
END_C_DECLS

#define C_NUMB_IIF_COMP 5

//////////////////////////////////////////////////////////////////////////////
/// Write NIT file of observer from luminance matrix of RGB components
/// @param nitfile - observer file(IN)
/// @param coldata - matrix of color luminance distribution(IN)
/// @param negvalue - how to process negative values
/// (-1 - keep as is, 0 - move to zero, 1 - reverse sign)(IN)
/// @return SUCCESS/FAILURE
OKAY WriteNITFile(const PathStr& nitfile,
    const TMatrix<Vect3f>& coldata, int negvalue)
{
    double white[XY], red[XY], green[XY], blue[XY];
    Str buf, buf1, buf2, layers, format, type;
    int i, j, k;
    float* pf, ** table;
    IIF* iif_file;
    INT64 rays;

    if ((iif_file = iif_open(nitfile.XData(), "w")) == NULL)
    {
        printf("\nIt is impossible to create observer file - %s",
            nitfile.XData());
        return FAILURE;
    }

    // Initialize resolution and component types
    layers = "lum red,lum gre,lum blu,lum acc,lum ray";
    format = "fffff";
    type = "LUMINANCE";
    if (iif_init_file(iif_file, coldata.NColumns(), coldata.NRows(),
        layers.XData(), format.XData()) != IIF_OK)
    {
        printf("\nIt is impossible to create observer file - %s",
            nitfile.XData());
        iif_close(iif_file);
        return FAILURE;
    }

    if (iif_put_var(iif_file, "FILE_TYPE", type.XData()) != IIF_OK)
    {
        printf("\nIt is impossible to create observer file - %s",
            nitfile.XData());
        iif_close(iif_file);
        return FAILURE;
    }

    if (iif_put_var(iif_file, "image pixel step", "1 1") != IIF_OK)
    {
        printf("\nIt is impossible to create observer file - %s",
            nitfile.XData());
        iif_close(iif_file);
        return FAILURE;
    }

    buf.Printf("1 1");
    if (iif_put_var(iif_file, "step size [m]", buf.XData()) != IIF_OK)
    {
        printf("\nIt is impossible to create observer file - %s",
            nitfile.XData());
        iif_close(iif_file);
        return FAILURE;
    }

    if (col_get_wrgb(white, red, green, blue) != COL_OK)
        ASSERT(FALSE)
    else
    {
        if (iif_put_var(iif_file, "GAMUT", "Yes") != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }

        buf.Printf(IF_V_FORMAT(IF_V_GAMUT), white[X]);
        if (iif_put_var(iif_file, "WHITE_X", buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }
        buf.Printf(IF_V_FORMAT(IF_V_GAMUT), white[Y]);
        if (iif_put_var(iif_file, "WHITE_Y", buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }

        buf.Printf(IF_V_FORMAT(IF_V_GAMUT), red[X]);
        if (iif_put_var(iif_file, "RED_X", buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }
        buf.Printf(IF_V_FORMAT(IF_V_GAMUT), red[Y]);
        if (iif_put_var(iif_file, "RED_Y", buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }

        buf.Printf(IF_V_FORMAT(IF_V_GAMUT), green[X]);
        if (iif_put_var(iif_file, "GREEN_X", buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }
        buf.Printf(IF_V_FORMAT(IF_V_GAMUT), green[Y]);
        if (iif_put_var(iif_file, "GREEN_Y", buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }

        buf.Printf(IF_V_FORMAT(IF_V_GAMUT), blue[X]);
        if (iif_put_var(iif_file, "BLUE_X", buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }
        buf.Printf(IF_V_FORMAT(IF_V_GAMUT), blue[Y]);
        if (iif_put_var(iif_file, "BLUE_Y", buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }
        rays = 0;
        buf.Printf(IF_V_FORMAT(IF_V_RAY_NUMBER), rays);
        if (iif_put_var(iif_file, IF_V_NAME(IF_V_RAY_NUMBER),
            buf.XData()) != IIF_OK)
        {
            printf("\nIt is impossible to create observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            return FAILURE;
        }
    }

    /* Set data pointers */
    table = new float* [C_NUMB_IIF_COMP];
    if (table == NULL)
    {
        printf("\nMemory allocation error - LUX file");
        iif_close(iif_file);
        delete[] table;
        return FAILURE;
    }
    table[R] = new float[coldata.NColumns() * C_NUMB_IIF_COMP];

    if (table[R] == NULL)
    {
        iif_close(iif_file);
        printf("\nMemory allocation error - NIT file");
        return FAILURE;
    }
    table[G] = table[R] + G * coldata.NColumns();
    table[B] = table[R] + B * coldata.NColumns();
    table[B + 1] = table[R] + (B + 1) * coldata.NColumns();
    table[B + 2] = table[R] + (B + 2) * coldata.NColumns();

    /* Print image */
    for (j = 0; j < coldata.NRows(); j++)
    {
        for (k = 0; k < RGB; k++)
        {
            pf = table[k];
            for (i = 0; i < coldata.NColumns(); i++)
            {
                if (negvalue > 0)
                    *pf++ = (float)Abs(coldata(j, i)[k]);
                else if (negvalue == 0 && coldata(j, i)[k] < 0)
                    *pf++ = (float)0;
                else
                    *pf++ = (float)(coldata(j, i)[k]);
            }
        }
        for (; k < C_NUMB_IIF_COMP; k++)
        {
            pf = table[k];
            for (i = 0; i < coldata.NColumns(); i++)
                *pf++ = (float)0;
        }

        if (iif_write_line(iif_file, (void*)table, j, -1, 0) != IIF_OK)
        {
            printf("\nIt is impossible to write observer file - %s",
                nitfile.XData());
            iif_close(iif_file);
            delete[] table[R];
            delete[] table;
            return FAILURE;
        }
    }

    // Close iif_file 
    delete[] table[R];
    delete[] table;
    iif_close(iif_file);
    return SUCCESS;
} // End of WriteNITFile()

//////////////////////////////////////////////////////////////////////////
/// Embree device error callback function
void DeviceErrorFunction(void* userPtr, enum RTCError error, const char* str)
{
    printf("Device error %d: %s\n", error, str);
}

//////////////////////////////////////////////////////////////////////////
/// Create a box object for Embree
/// @param[in] device Embree device
/// @param[in] p Box origin points
/// @param[in] size Box sizes along 3 axis
/// @param[in] tr Transformation matrix to apply to created geometry
/// @return Embree geometry object containing box
RTCGeometry CreateBox(RTCDevice device, Point3f p, Vect3f size, Matrix43f tr)
{
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vertices = (float*)rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT3,
        3 * sizeof(float),
        8);

    unsigned* indices = (unsigned*)rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_INDEX,
        0,
        RTC_FORMAT_UINT3,
        3 * sizeof(unsigned),
        12);
    if (vertices && indices)
    {
        Point3f points[8];
        points[0] = p;
        points[1] = p + Vect3f(size.x, 0, 0);
        points[2] = p + Vect3f(size.x, 0, size.z);
        points[3] = p + Vect3f(0, 0, size.z);
        points[4] = p + Vect3f(0, size.y, size.z);
        points[5] = p + Vect3f(0, size.y, 0);
        points[6] = p + Vect3f(size.x, size.y, 0);
        points[7] = p + size;

        for (int i = 0; i < 8; i++)
        {
            tr.PointTransform(points[i]);
            vertices[i * 3] = points[i].x;
            vertices[i * 3 + 1] = points[i].y;
            vertices[i * 3 + 2] = points[i].z;
        }

        int i = 0;
        indices[i++] = 0;
        indices[i++] = 1;
        indices[i++] = 2;

        indices[i++] = 2;
        indices[i++] = 3;
        indices[i++] = 0;

        indices[i++] = 4;
        indices[i++] = 5;
        indices[i++] = 0;

        indices[i++] = 0;
        indices[i++] = 3;
        indices[i++] = 4;

        indices[i++] = 7;
        indices[i++] = 6;
        indices[i++] = 5;

        indices[i++] = 5;
        indices[i++] = 4;
        indices[i++] = 7;

        indices[i++] = 3;
        indices[i++] = 2;
        indices[i++] = 7;

        indices[i++] = 7;
        indices[i++] = 4;
        indices[i++] = 3;

        indices[i++] = 7;
        indices[i++] = 2;
        indices[i++] = 1;//

        indices[i++] = 1;
        indices[i++] = 6;
        indices[i++] = 7;//

        indices[i++] = 0;
        indices[i++] = 5;
        indices[i++] = 6;

        indices[i++] = 6;
        indices[i++] = 1;
        indices[i++] = 0;
    }

    rtcCommitGeometry(geom);

    return geom;
}

/// Which box faces to omit flags
enum OmitFace
{
    /// Omit nothing
    OMIT_NONE = 0,
    /// Omit +X face
    OMIT_X_POS = 1,
    // Omit -X face
    OMIT_X_NEG = 2,
    /// Omit +Y face
    OMIT_Y_POS = 4,
    /// Omit -Y face
    OMIT_Y_NEG = 8,
    /// Omit +Z face
    OMIT_Z_POS = 16,
    /// Omit -Z face
    OMIT_Z_NEG = 32
};

//////////////////////////////////////////////////////////////////////////
/// Create a box object for Embree. Some faces may be omitted.
/// @param[in] device Embree device
/// @param[in] p Box origin points
/// @param[in] size Box sizes along 3 axis
/// @param[in] omit A mask of which faces to omit
/// @param[in] tr Transformation matrix to apply to created geometry
/// @return Embree geometry object containing box
RTCGeometry CreateBoxOmit(RTCDevice device, Point3f p, Vect3f size, int omit, Matrix43f tr)
{
    int omits[6] = { 1, 2, 4, 8, 16, 32 };

    int omit_count = 0;
    for (int i = 0; i < 6; i++)
        if (omit & omits[i])
            omit_count++;

    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vertices = (float*)rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT3,
        3 * sizeof(float),
        8);

    unsigned* indices = (unsigned*)rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_INDEX,
        0,
        RTC_FORMAT_UINT3,
        3 * sizeof(unsigned),
        12 - omit_count * 2);
    if (vertices && indices)
    {
        Point3f points[8];
        points[0] = p;
        points[1] = p + Vect3f(size.x, 0, 0);
        points[2] = p + Vect3f(size.x, 0, size.z);
        points[3] = p + Vect3f(0, 0, size.z);
        points[4] = p + Vect3f(0, size.y, size.z);
        points[5] = p + Vect3f(0, size.y, 0);
        points[6] = p + Vect3f(size.x, size.y, 0);
        points[7] = p + size;

        for (int i = 0; i < 8; i++)
        {
            tr.PointTransform(points[i]);
            vertices[i * 3] = points[i].x;
            vertices[i * 3 + 1] = points[i].y;
            vertices[i * 3 + 2] = points[i].z;
        }

        int i = 0;
        if (!(omit & OMIT_Y_NEG))
        {
            indices[i++] = 0;
            indices[i++] = 1;
            indices[i++] = 2;

            indices[i++] = 2;
            indices[i++] = 3;
            indices[i++] = 0;
        }

        if (!(omit & OMIT_X_NEG))
        {
            indices[i++] = 4;
            indices[i++] = 5;
            indices[i++] = 0;

            indices[i++] = 0;
            indices[i++] = 3;
            indices[i++] = 4;
        }

        if (!(omit & OMIT_Y_POS))
        {
            indices[i++] = 7;
            indices[i++] = 6;
            indices[i++] = 5;

            indices[i++] = 5;
            indices[i++] = 4;
            indices[i++] = 7;
        }

        if (!(omit & OMIT_Z_POS))
        {
            indices[i++] = 3;
            indices[i++] = 2;
            indices[i++] = 7;

            indices[i++] = 7;
            indices[i++] = 4;
            indices[i++] = 3;
        }

        if (!(omit & OMIT_X_POS))
        {
            indices[i++] = 7;
            indices[i++] = 2;
            indices[i++] = 1;

            indices[i++] = 1;
            indices[i++] = 6;
            indices[i++] = 7;
        }

        if (!(omit & OMIT_Z_NEG))
        {
            indices[i++] = 0;
            indices[i++] = 5;
            indices[i++] = 6;

            indices[i++] = 6;
            indices[i++] = 1;
            indices[i++] = 1;
        }
    }

    rtcCommitGeometry(geom);

    return geom;
}

//////////////////////////////////////////////////////////////////////////
/// Subdivide sphere face
/// @param[in] v1 First vertex of face to subdivide
/// @param[in] v2 Second vertex of face to subdivide
/// @param[in] v3 Third vertex of face to subdivide
/// @param[in] i1 Index of the first vertex
/// @param[in] i2 Index of the second vertex
/// @param[in] i3 Index of the third vertex
/// @param[in, out] sphere_pints Array of vertexes
/// @param[in, out] sphere_indices Array of triangle indices
/// @param[in] depth Subdivision depth
void SphereSubdivide(Vect3f v1, Vect3f v2, Vect3f v3, const unsigned i1, const unsigned i2, const unsigned i3,
    TArray<Vect3f>& sphere_points, TArray<Vect3u>& sphere_indices, unsigned int depth) {
    if (depth == 0)
    {
        sphere_indices.Add(Vect3u(i1, i2, i3));
        return;
    }

    Vect3f v12 = (v1 + v2).Normalize();
    Vect3f v23 = (v2 + v3).Normalize();
    Vect3f v31 = (v3 + v1).Normalize();

    unsigned i12 = sphere_points.Length();
    unsigned i23 = i12 + 1;
    unsigned i31 = i23 + 1;

    sphere_points.Add(v12);
    sphere_points.Add(v23);
    sphere_points.Add(v31);

    SphereSubdivide(v1, v12, v31, i1, i12, i31, sphere_points, sphere_indices, depth - 1);
    SphereSubdivide(v2, v23, v12, i2, i23, i12, sphere_points, sphere_indices, depth - 1);
    SphereSubdivide(v3, v31, v23, i3, i31, i23, sphere_points, sphere_indices, depth - 1);
    SphereSubdivide(v12, v23, v31, i12, i23, i31, sphere_points, sphere_indices, depth - 1);
}

//////////////////////////////////////////////////////////////////////////
/// Start subdividing a sphere
/// @param[out] sphere_pints Array of vertexes
/// @param[out] sphere_indices Array of triangle indices
/// @param[in] depth Subdivision depth
void InitSphere(TArray<Vect3f>& sphere_points, TArray<Vect3u>& sphere_indices, unsigned int depth)
{
    const float X0 = 0.525731112119133606f;
    const float Z0 = 0.850650808352039932f;
    sphere_points.Truncate();
    sphere_points.Add({ -X0, 0.0,  Z0 });
    sphere_points.Add({ X0, 0.0,  Z0 });
    sphere_points.Add({ -X0, 0.0, -Z0 });
    sphere_points.Add({ X0, 0.0, -Z0 });
    sphere_points.Add({ 0.0,  Z0,  X0 });
    sphere_points.Add({ 0.0,  Z0, -X0 });
    sphere_points.Add({ 0.0, -Z0,  X0 });
    sphere_points.Add({ 0.0, -Z0, -X0 });
    sphere_points.Add({ Z0,  X0, 0.0 });
    sphere_points.Add({ -Z0,  X0, 0.0 });
    sphere_points.Add({ Z0, -X0, 0.0 });
    sphere_points.Add({ -Z0, -X0, 0.0 });


    int tindices[20][3] =
    {
      {0, 4, 1},    { 0, 9, 4 },  { 9, 5, 4 },  { 4, 5, 8 },  { 4, 8, 1 },
      { 8, 10, 1 }, { 8, 3, 10 }, { 5, 3, 8 },  { 5, 2, 3 },  { 2, 7, 3 },
      { 7, 10, 3 }, { 7, 6, 10 }, { 7, 11, 6 }, { 11, 0, 6 }, { 0, 1, 6 },
      { 6, 1, 10 }, { 9, 0, 11 }, { 9, 11, 2 }, { 9, 2, 5 },  { 7, 2, 11 }
    };
    for (int i = 0; i < 20; i++)
        SphereSubdivide(sphere_points[tindices[i][0]], sphere_points[tindices[i][1]], sphere_points[tindices[i][2]],
            tindices[i][0], tindices[i][1], tindices[i][2],
            sphere_points, sphere_indices, depth);
}

//////////////////////////////////////////////////////////////////////////
/// Create a sphere object for Embree
/// @param[in] device Embree device
/// @param[in] center Sphere center
/// @param[in] radius Sphere radius
/// @param[in] tr Transformation matrix to apply to created geometry
/// @return Embree geometry object containing sphere
RTCGeometry CreateSphere(RTCDevice device, Point3f center, float radius, unsigned int depth, Matrix43f tr)
{
    TArray<Vect3f> sphere_normals;
    TArray<Vect3u> sphere_indices;

    InitSphere(sphere_normals, sphere_indices, depth);

    TArray<Point3f> sphere_points;
    for (int i = 0; i < sphere_normals.Length(); i++)
    {
        sphere_points.Add(center + sphere_normals[i] * radius);
        tr.PointTransform(sphere_points[i]);
    }

    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vertices = (float*)rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT3,
        3 * sizeof(float),
        sphere_points.Length());

    unsigned* indices = (unsigned*)rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_INDEX,
        0,
        RTC_FORMAT_UINT3,
        3 * sizeof(unsigned),
        sphere_indices.Length());

    if (vertices && indices)
    {
        for (int i = 0; i < sphere_points.Length(); i++)
        {
            vertices[i * 3] = sphere_points[i].x;
            vertices[i * 3 + 1] = sphere_points[i].y;
            vertices[i * 3 + 2] = sphere_points[i].z;
        }

        for (int i = 0; i < sphere_indices.Length(); i++)
        {
            indices[i * 3] = sphere_indices[i].x;
            indices[i * 3 + 1] = sphere_indices[i].y;
            indices[i * 3 + 2] = sphere_indices[i].z;
        }
    }

    rtcCommitGeometry(geom);

    return geom;
}

// ������� �������

float func1(float x) {
    return x * x;
}

float func2(float x) {
    return x * x;
}

float monteCarlo1(float A, float C, int N) {
    float sum = 0;
    Rnd random = Rnd();

    for (int i = 0; i < N; i++) {
        float r = random.DRnd() * (C - A) + A;

        sum += func1(r);
    }

    sum *= (C - A);
    sum /= N;

    return sum;
}

float monteCarlo2(float A, float C, int N, int k) {
    float sum2 = 0;
    float sum1 = 0;
    Rnd random = Rnd();
    int c = N / k;

    for (int j = 0; j < k; j++) {
        for (int i = 0; i < c; i++) {
            float r = random.DRnd() * (C - A) / k + A + j * ((C - A) / k);

            sum1 += func1(r);
        }

        sum1 *= (C - A);
        sum1 /= c;
        sum2 += sum1;
        sum1 = 0;
    }

    return sum2 / k;
}

float monteCarlo3(float A, float C, int N, int L) {
    float sum = 0;
    float sq = (C * C * C - A * A * A) / 3;;
    Rnd random = Rnd();
    float Ay = func1(A);
    float Cy = func1(C);
    for (int i = 0; i < L; i++) {
        float y = random.DRnd() * (Cy - Ay) + Ay;
        float x = sqrt(y);
        sum += func1(x) / func2(x) * sq;
    }
    return sum / L;
}

class datchick {
    double x, y, z, r;

public:
    datchick(double a, double b, double c, double d) {
        x = a;
        y = b;
        z = c;
        r = d;
    }

    int popal(double a, double b, double c) {
        if ((a - x < r) && (a - x > -r))
            if ((b - y < r) && (b - y > -r))
                if ((c - z < r) && (c - z > -r))
                    if ((a - x) * (a - x) + (b - y) * (b - y) + (c - z) * (c - z) <= r * r)
                        return 1;
        return 0;
    }
};

void RectRaspledelenie() { // rectangle
    int N = 1000000;
    double a = 5, b = 3;
    datchick d1 = datchick(1, 1, 0, 1);
    datchick d2 = datchick(3, 2, 0, 1);

    Rnd random = Rnd();

    int dat1 = 0, dat2 = 0;

    for (int i = 0; i < N; i++) {
        double ar = random.DRnd() * a;
        double br = random.DRnd() * b;

        if (d1.popal(ar, br, 0))
            dat1++;
        if (d2.popal(ar, br, 0))
            dat2++;
    }

    double res = N * 3.14 / (a * b);

    printf("Rectangle\nreal result: %f \ndatchick 1: %d\ndatchick 2: %d\n\n", res, dat1, dat2);
}

void TriangleRaspledelenie() { // triangle
    int N = 1000000;
    double ax = 5, ay = 7, bx = 3, by = -3;
    datchick d1 = datchick(2, 0, 0, 1);
    datchick d2 = datchick(2, 1, 0, 1);

    Rnd random = Rnd();

    int dat1 = 0, dat2 = 0;

    for (int i = 0; i < N; i++) {
        double ir = random.DRnd();
        double jr = random.DRnd();

        if (ir + jr > 1) {
            ir = 1 - ir;
            jr = 1 - jr;
        }

        if (d1.popal(ax * ir + bx * jr, ay * ir + by * jr, 0))
            dat1++;
        if (d2.popal(ax * ir + bx * jr, ay * ir + by * jr, 0))
            dat2++;
    }

    double res = N * 3.14 / (0.5 * Abs((ax * by - bx * ay)));

    printf("Triangle\nreal result: %f \ndatchick 1: %d\ndatchick 2: %d\n\n", res, dat1, dat2);
}

void CircleRaspledelenie() { // circle
    int N = 1000000;
    double r = 5;
    datchick d1 = datchick(0, 0, 0, 1);
    datchick d2 = datchick(2, 2, 0, 1);
    datchick d3 = datchick(-3, -2, 0, 1);

    Rnd random = Rnd();

    int dat1 = 0, dat2 = 0, dat3 = 0;

    for (int i = 0; i < N; i++) {
        double xr = random.DRnd() * 2 * r - r;
        double yr = random.DRnd() * 2 * r - r;

        if (xr * xr + yr * yr > r * r) {
            i--;
            continue;
        }

        if (d1.popal(xr, yr, 0))
            dat1++;
        if (d2.popal(xr, yr, 0))
            dat2++;
        if (d3.popal(xr, yr, 0))
            dat3++;
    }

    double res = N / r / r;

    printf("Circle\nreal result: %f \ndatchick 1: %d\ndatchick 2: %d\ndatchick 3: %d\n\n", res, dat1, dat2, dat3);
}

void PoverhnoctSphereRaspledelenie() { // sphere
    int N = 1000000;
    double r = 5;

    double d1h1 = 1, d1h2 = 2;
    double d1phi1 = 0, d1phi2 = PI / 2;
    double d2h1 = -3, d2h2 = -2;
    double d2phi1 = PI, d2phi2 = PI * 3 / 2;

    Rnd random = Rnd();

    int dat1 = 0, dat2 = 0;

    for (int i = 0; i < N; i++) {
        double hr = random.DRnd() * 2 * r - r;
        double phir = random.DRnd() * 2 * PI;

        if ((hr > d1h1 && hr < d1h2) && (phir > d1phi1 && phir < d1phi2))
            dat1++;
        if ((hr > d2h1 && hr < d2h2) && (phir > d2phi1 && phir < d2phi2))
            dat2++;
    }

    double res = N * (2 * PI * r) / 4 / (4 * PI * r * r);

    printf("Poverchnost sphere\nreal result: %f \ndatchick 1: %d\ndatchick 2: %d\n\n", res, dat1, dat2);
}

void VnutriSphereRaspledelenie() { // circle
    int N = 1000000;
    double r = 5;
    datchick d1 = datchick(0, 0, 0, 1);
    datchick d2 = datchick(2, 2, 2, 1);
    datchick d3 = datchick(-3, -2, -1, 1);

    Rnd random = Rnd();

    int dat1 = 0, dat2 = 0, dat3 = 0;

    for (int i = 0; i < N; i++) {
        double xr = random.DRnd() * 2 * r - r;
        double yr = random.DRnd() * 2 * r - r;
        double zr = random.DRnd() * 2 * r - r;

        if (xr * xr + yr * yr + zr * zr > r * r) {
            i--;
            continue;
        }

        if (d1.popal(xr, yr, zr))
            dat1++;
        if (d2.popal(xr, yr, zr))
            dat2++;
        if (d3.popal(xr, yr, zr))
            dat3++;
    }

    double res = N / r / r / r;

    printf("Vnutri sphere\nreal result: %f \ndatchick 1: %d\ndatchick 2: %d\ndatchick 3: %d\n\n", res, dat1, dat2, dat3);
}

//classes: light, camera(+�����), ray, opt prop, geometry -> scene -> render

Rnd random = Rnd();

struct vf {
    size_t v;
    size_t vn;
    size_t f;
};

vf count(FILE* file) {
    char c = fgetc(file);
    size_t v = 0, vn = 0, f = 0;

    while (c != EOF) {
        switch (c)
        {
        case '#': {
            do {
                c = fgetc(file);
            } while (c != '\n');
            c = fgetc(file);
            continue;
        }
        case 'v': {
            c = fgetc(file);
            if (c == ' ')
                v++;
            if (c == 'n')
                vn++;
            /*else {
                printf("%c", c);
            }*/
            do {
                c = fgetc(file);
            } while (c != '\n');
            c = fgetc(file);
            continue;
        }
        case 'f': {
            f++;
            printf("%zu\n", f);
            do {
                c = fgetc(file);
            } while (c != '\n');
            c = fgetc(file);
            continue;
        }
        default: {
            do {
                c = fgetc(file);
            } while (c != '\n');
            c = fgetc(file);
        }
        }
    }
    //printf("%d, %d", v, f);
    return { v, vn, f };
}

int read_int(FILE* f) {
    char c = fgetc(f);
    int res = 0;
    int slash = 0;

    while (c != ' ' && c != '\n' && c != EOF) {
        if (slash) {
            c = fgetc(f);
            continue;
        }
        if (c == '/') {
            slash = 1;
            continue;
        }
        switch (c) {
        case '0': {
            res *= 10;
            res += 0;
            c = fgetc(f);
            break;
        }
        case '1': {
            res *= 10;
            res += 1;
            c = fgetc(f);
            break;
        }
        case '2': {
            res *= 10;
            res += 2;
            c = fgetc(f);
            break;
        }
        case '3': {
            res *= 10;
            res += 3;
            c = fgetc(f);
            break;
        }
        case '4': {
            res *= 10;
            res += 4;
            c = fgetc(f);
            break;
        }
        case '5': {
            res *= 10;
            res += 5;
            c = fgetc(f);
            break;
        }
        case '6': {
            res *= 10;
            res += 6;
            c = fgetc(f);
            break;
        }
        case '7': {
            res *= 10;
            res += 7;
            c = fgetc(f);
            break;
        }
        case '8': {
            res *= 10;
            res += 8;
            c = fgetc(f);
            break;
        }
        case '9': {
            res *= 10;
            res += 9;
            c = fgetc(f);
            break;
        }
        }
    }
    return res;
}

float read_float(FILE* f) {
    char c = fgetc(f);
    float res = 0;
    int tens = 0;
    int minus = 0;

    while (c == ' ') {
        c = fgetc(f);
    }

    if (c == '-') {
        minus = 1;
        c = fgetc(f);
    }

    while (c != ' ' && c != '\n' && c != EOF) {
        switch (c) {
        case '0': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 0;
            c = fgetc(f);
            break;
        }
        case '1': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 1;
            c = fgetc(f);
            break;
        }
        case '2': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 2;
            c = fgetc(f);
            break;
        }
        case '3': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 3;
            c = fgetc(f);
            break;
        }
        case '4': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 4;
            c = fgetc(f);
            break;
        }
        case '5': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 5;
            c = fgetc(f);
            break;
        }
        case '6': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 6;
            c = fgetc(f);
            break;
        }
        case '7': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 7;
            c = fgetc(f);
            break;
        }
        case '8': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 8;
            c = fgetc(f);
            break;
        }
        case '9': {
            if (tens)
                tens *= 10;
            res *= 10;
            res += 9;
            c = fgetc(f);
            break;
        }
        case '.': {
            tens = 1;
            c = fgetc(f);
            break;
        }
        }
    }
    if (minus)
        return -res / tens;
    return res / tens;
}

Vect3f* VN;
Vect3i* normals;

RTCGeometry ReadOBJ(RTCDevice device, char* objfile, float x)
{
    FILE* f = fopen(objfile, "r");
    vf tmp = count(f);
    printf("%zu %zu %zu\n", tmp.f, tmp.v, tmp.vn);
    fclose(f);
    f = fopen(objfile, "r");

    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vertices = (float*)rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT3,
        3 * sizeof(float),
        tmp.v);

    /*float* normals = (float*)rtcSetNewGeometryBuffer(geom, 
        RTC_BUFFER_TYPE_NORMAL,
        0,
        RTC_FORMAT_FLOAT3,
        3 * sizeof(float), 
        tmp.vn);*/

    VN = new Vect3f[tmp.vn];

    normals = new Vect3i[tmp.f * 2];

    unsigned* indices = (unsigned*)rtcSetNewGeometryBuffer(geom,
        RTC_BUFFER_TYPE_INDEX,
        0,
        RTC_FORMAT_UINT3,
        3 * sizeof(unsigned),
        tmp.f * 2);

    if (vertices && indices)
    {
        char ch = fgetc(f);
        size_t i = 0;
        size_t in = 0;
        size_t j = 0;

        while (ch != EOF) {
            switch (ch) {
            case 'v': {
                ch = fgetc(f);
                if (ch == ' ') {
                    float a = read_float(f), b = read_float(f), c = read_float(f);
                    vertices[i * 3] = a;
                    vertices[i * 3 + 1] = b;
                    vertices[i * 3 + 2] = c - x;
                    if (i < 10)
                        printf("%f %f %f\n", a, b, c);
                    i++;
                }
                else if (ch == 'n') {
                    float a = read_float(f), b = read_float(f), c = read_float(f);
                    VN[in].x = a;
                    VN[in].y = b;
                    VN[in].z = c;
                    if (in < 10)
                        printf("%f %f %f\n", a, b, c);
                    VN[in].Normalize();
                    in++;
                }
                else {
                    while (ch != '\n') {
                        ch = fgetc(f);
                    }
                    ch = fgetc(f);
                }

                break;
            }
            case 'f': {
                ch = fgetc(f);
                if (ch != ' ')
                    break;
                int a = read_int(f), b = read_int(f), c = read_int(f), d = read_int(f);
                indices[j * 3] = a - 1;
                indices[j * 3 + 1] = b - 1;
                indices[j * 3 + 2] = c - 1;
                //printf("%d %d %d\n", a, b, c);
                normals[j].x = a - 1;
                normals[j].y = b - 1;
                normals[j].z = c - 1;
                j++;

                if (d != 0) {
                    indices[j * 3] = c - 1;
                    indices[j * 3 + 1] = d - 1;
                    indices[j * 3 + 2] = a - 1;

                    normals[j].x = c - 1;
                    normals[j].y = d - 1;
                    normals[j].z = a - 1;

                    j++;
                }
                else {
                }

                //printf("%d %d %d %d\n", a, b, c, d);
                break;
            }
            default: {
                //printf("hi\n");
                while (ch != '\n') {
                    ch = fgetc(f);
                }
            }
            }
            ch = fgetc(f);
        }
    }
    fclose(f);
    rtcCommitGeometry(geom);

    return geom;
}

class ray { // ����� ����
public:
    Point3f start; // ����� ������
    Vect3f dir; // �����������
    float tfar; // ����� ����
    Vect3f rgb; // ���� ����
    float n; // ���������� ����������� ����

public:
    ray() {
        n = 1;
    }
    ray(Point3f st, Point3f fin) { // �������� ���� �� ������ ������ � �����
        start = st;
        dir = Vect3f(fin.x - st.x, fin.y - st.y, fin.z - st.z);
        tfar = dir.Length();
        if (tfar != 0)
            dir.Normalize();
        rgb = Vect3f(1, 1, 1);
        rgb.Normalize();
        n = 1.0;
    }
    ray(Point3f st, Vect3f direction) { // �������� ���� �� ������ ������ � �����
        start = st;
        dir = direction;
        tfar = dir.Length();
        if (tfar != 0)
            dir.Normalize();
        rgb = Vect3f(1, 1, 1);
        rgb.Normalize();
        n = 1.0;
    }
    Vect3f mirror_color(Vect3f E) { // ���� �������
        Vect3f res = Vect3f(rgb.x * E.x, rgb.y * E.y, rgb.z * E.z);
        return res;
    }
};

//class light {
//    Vect3f Intensity;
//
//public:
//    Vect3f illumination(Point3f p, Vect3f N) {
//        return Vect3f(0);
//    }
//};

//class light{ // ����� ��������� �����
//public:
//    Point3f position; // ��������� 
//    Vect3f Intensity; // ������������� �� ���� ������� ����� ������� 
//
//public:
//    light() {
//        position = Point3f(0, 0, 0);
//        Intensity = Vect3f(1000, 1000, 1000);
//    }
//    light(float x, float y, float z, float r, float g, float b) {
//        position = Point3f(x, y, z);
//        Intensity = Vect3f(r, g, b);
//    }
//    Vect3f illumination(Point3f p, Vect3f N) { // ������� ������������
//        ray r = ray(p, position); //��� �� ����� �� ������� �� ��������� �����
//        float d = r.tfar; //���������� �� ��������� �����
//        float cos = r.dir.x * N.x + r.dir.y * N.y + r.dir.z * N.z; // ������� � ������ ���� ���������
//        if (d == 0)
//            return 0;
//        else {
//            /*if (cos < 0)
//                return Vect3f(0);*/
//            if (cos < 0)
//                cos = -cos;
//            Vect3f res = Intensity * cos / d / d;
//            return res;
//        }
//    }
//};

Vect3f vectMult(Vect3f a, Vect3f b) {
    float x = a.y * b.z - a.z * b.y;
    float y = a.x * b.z - a.z * b.x;
    float z = a.x * b.y - a.y * b.x;

    return Vect3f(x, y, z);
}

float scalarMult(Vect3f a, Vect3f b) {
    float x = a.x * b.x;
    float y = a.y * b.y;
    float z = a.z * b.z;

    return x + y + z;
}

class light { // ����� ��������� �����
public:
    //Point3f position; // ��������� 
    //Vect3f Intensity; // ������������� �� ���� ������� ����� ������� 
    Point3f p1, p2, p3;
    Vect3f Lightness;// �������
    Vect3f Normal;
    float S;

public:
    light() {
        p1 = Point3f(0, 0, 0);
        p2 = Point3f(0, 0, 0);
        p3 = Point3f(0, 0, 0);
        Lightness = Vect3f(1000, 1000, 1000);
        Normal = Vect3f(1, 0, 0);
        S = 1;
    }
    light(Point3f x, Point3f y, Point3f z, float red, float green, float blue) {
        p1 = x;
        p2 = y;
        p3 = z;
        Lightness = Vect3f(red, green, blue);
        //Normal = Vect3f(p2 - p1).operator*(Vect3f(p3 - p1));
        Normal = vectMult(Vect3f(p2 - p1), Vect3f(p3 - p1));
        Normal.Normalize();
        float p = (Vect3f(p2 - p1).Length() + Vect3f(p3 - p2).Length() + Vect3f(p1 - p3).Length()) / 2;
        float a = Vect3f(p2 - p1).Length();
        float b = Vect3f(p3 - p2).Length();
        float c = Vect3f(p1 - p3).Length();
        S = Sqrt(p * (p - a) * (p - b) * (p - c));
    }

    Point3f randomP() {
        float ir = random.DRnd();
        float jr = random.DRnd();

        if (ir + jr > 1) {
            ir = 1 - ir;
            jr = 1 - jr;
        }

        Vect3f v1 = ir * (p2 - p1);
        Vect3f v2 = jr * (p3 - p1);
        Vect3f res = v1 + v2;

        return p1 + res;
    }

    RTCGeometry lightGeom(RTCDevice device)
    {
        RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
        float* vertices = (float*)rtcSetNewGeometryBuffer(geom,
            RTC_BUFFER_TYPE_VERTEX,
            0,
            RTC_FORMAT_FLOAT3,
            3 * sizeof(float),
            3);

        unsigned* indices = (unsigned*)rtcSetNewGeometryBuffer(geom,
            RTC_BUFFER_TYPE_INDEX,
            0,
            RTC_FORMAT_UINT3,
            3 * sizeof(unsigned),
            1);

        if (vertices && indices)
        {
            vertices[0] = p1.x;
            vertices[1] = p1.y;
            vertices[2] = p1.z;

            vertices[3] = p2.x;
            vertices[4] = p2.y;
            vertices[5] = p2.z;

            vertices[6] = p3.x;
            vertices[7] = p3.y;
            vertices[8] = p3.z;

            indices[0] = 0;
            indices[1] = 1;
            indices[2] = 2;
        }

        rtcCommitGeometry(geom);

        return geom;
    }

    Vect3f illumination(Point3f p, Point3f g, Vect3f N) { // ������� ������������
        ray r = ray(p, g); //��� �� ����� �� ������� �� ��������� �����
        float d = r.tfar; //���������� �� ��������� �����
        float cos_gamma = r.dir.x * N.x + r.dir.y * N.y + r.dir.z * N.z; // ������� � ������ ���� ���������
        float cos_tetha = r.dir.x * Normal.x + r.dir.y * Normal.y + r.dir.z * Normal.z; // ������� � ������ ���� ���������
        if (d == 0)
            return 0;
        else {
            /*if (cos < 0)
                return Vect3f(0);*/
            if (cos_gamma < 0)
                cos_gamma = -cos_gamma;

            if (cos_tetha < 0)
                cos_tetha = -cos_tetha;

            Vect3f res = Lightness * S * cos_tetha * cos_gamma / d / d;
            //Vect3f res = Lightness * S * cos_gamma / d / d;
            return res;
        }
    }

    float F() {
        float d = Lightness.Length();
        return d * S * PI;
    }
};

ray RandomRay(Point3f a, Vect3f n) { // sphere
    double hr = random.DRnd() * 2 * -1;
    double phir = random.DRnd() * 2 * PI;

    Vect3f r;

    if (n.z != 0) {
        r.x = n.y;
        r.y = -n.x;
        r.z = 0;
    }
    else {
        r.x = 0;
        r.y = -n.z;
        r.z = n.y;
    }

    Matrix43f tr = Matrix43f(1, 1, 1);
    tr.RotationAxis(n, phir);
    tr.VectorTransform(r);

    return ray(a, n * hr + r);
}

class appearance {
public:
    appearance(){}
    virtual ~appearance(){}

    virtual double choose_event(double r, ray* ray) const = 0;

    virtual int diffusion(ray* ray, Vect3f normal) const = 0;

    virtual Vect3f luminance() const = 0;
};

class mirror_reflection : public appearance {
    Vect3f ks;
public:
    mirror_reflection() : appearance(){
        ks = 0;
    }

    mirror_reflection(Vect3f x) : appearance(){ //����� ������
        ks = x;
    }
    //����������� �����������
    //�������� �����

    virtual ~mirror_reflection(){}

    virtual double choose_event(double r, ray* ray) const
    {
        double x = DotProd((*ray).rgb, ks);
        r -= x;
        return r;
    }

    virtual int diffusion(ray* ray, Vect3f normal) const { //
        (*ray).dir = (*ray).dir - 2 * scalarMult((*ray).dir, normal) * normal;
        //printf("hi1\n");
        return 0; //����� enum        
    }

    virtual Vect3f luminance() const {
        return Vect3f(0);
    };
};

class diffuse_reflection : public appearance {
    Vect3f kd;
public:
    diffuse_reflection() : appearance() {
        kd = 0;
    }

    diffuse_reflection(Vect3f x) : appearance() { //����� ������
        kd = x;
    }
    //����������� �����������
    //�������� �����

    virtual ~diffuse_reflection() {}

    double choose_event(double r, ray* ray) const{
        double x = DotProd((*ray).rgb, kd);
        r -= x;
        return r;
    }
    virtual int diffusion(ray* ray, Vect3f normal) const {
        //printf("hi2\n");
        (*ray).dir = RandomRay(Point3f(0), normal).dir;
        return 1;
    }

    virtual Vect3f luminance() const {
        return Vect3f(0);
    };
};

class mirror_refraction : public appearance {
    Vect3f kts;
    float n1, n2;
public:
    mirror_refraction() : appearance() {
        kts = 0;
    }

    mirror_refraction(Vect3f x, float _n1, float _n2) : appearance() { //����� ������
        kts = x;
        n1 = _n1;
        n2 = _n2;
    }
    //����������� �����������
    //�������� �����

    virtual ~mirror_refraction() {}

    double choose_event(double r, ray* ray) const {
        double x = DotProd((*ray).rgb, kts);
        r -= x;
        return r;
    }
    virtual int diffusion(ray* ray, Vect3f normal) const {
        /*float refr_ind;
        if ((*ray).n == n1) {
            (*ray).n = n2;
            refr_ind = n1 / n2;
        }            
        else {
            (*ray).n = n1;
            refr_ind = n2 / n1;
        }

        float cos_in = Abs(scalarMult((*ray).dir, normal));
        float sin_in = 1 - cos_in * cos_in;
        sin_in = sin_in < 0 ? 0 : Sqrt(sin_in);

        float sin_out = sin_in / refr_ind;
        float cos_out = 1 - sin_out * sin_out;

        if (cos_out <= 0) {
            (*ray).dir = (*ray).dir - 2 * scalarMult((*ray).dir, normal) * normal;
            return 2;
        }

        cos_out = Sqrt(cos_out);

        float weight = cos_out * refr_ind - cos_in;

        if (scalarMult((*ray).dir, normal) < 0)
            (*ray).dir -= normal * weight;
        else
            (*ray).dir += normal * weight;

        (*ray).dir.Normalize();*/

        float eta; // eta = in_IOR/out_IOR

        if ((*ray).n == n1) {
            (*ray).n = n2;
            eta = n1 / n2;
        }
        else {
            (*ray).n = n1;
            eta = n2 / n1;
        }

        float cos_theta = scalarMult(normal, (*ray).dir);

        if (cos_theta < 0)
        {
            cos_theta *= -1.0f;
            normal *= -1.0f;
            eta = 1.0f / eta;

        }

        float k = 1.0f - eta * eta * (1.0 - cos_theta * cos_theta);

        if (k >= 0.0f)
            (*ray).dir = eta * (*ray).dir + (eta * cos_theta - sqrt(k)) * normal;

        (*ray).dir.Normalize();

        return 2; //����� enum
    }

    virtual Vect3f luminance() const {
        return Vect3f(0);
    };
};

class diffuse_refraction : public appearance {
    Vect3f ktd;
public:
    diffuse_refraction() : appearance() {
        ktd = 0;
    }

    diffuse_refraction(Vect3f x, int n1, int n2) : appearance() { //����� ������
        ktd = x;
    }
    //����������� �����������
    //�������� �����

    virtual ~diffuse_refraction() {}

    double choose_event(double r, ray* ray) const {
        double x = DotProd((*ray).rgb, ktd);
        r -= x;
        return r;
    }
    virtual int diffusion(ray* ray, Vect3f normal) const {
        (*ray).dir = - RandomRay(Point3f(0), normal).dir;
        return 3;
    }

    virtual Vect3f luminance() const {
        return Vect3f(0);
    };
};

class optprop { // ����� ����������� ��������
public:
    Vect3f rgb; // ���� �����������
    //int is_mirror; // �������� �� ����������� ��������
    int is_light; // �������� �� ����������� ���������� �����
    //float ks, kd, kts, ktd; //������ ���� �������
    float n;
    int count_event = 4; //���������� �������-�����������
    appearance *app[4]; //������� ����������

public:
    optprop() {
        rgb = Vect3f(255, 255, 255);
        rgb = rgb.Normalize(); // ����������� �� ������ �������
        is_light = 0;
        count_event = 4;
        app[0] = new mirror_reflection();
        app[1] = new diffuse_reflection();
        app[2] = new mirror_refraction();
        app[3] = new diffuse_refraction();
    }

    optprop(float r, float g, float b, int light, float _ks, float _kd, float _kts, float _ktd, float _n1, float _n2) {
        rgb = Vect3f(r, g, b);
        rgb = rgb.Length() == 0 ? rgb : rgb.Normalize(); // ����������� �� ������ �������
        is_light = light;
        count_event = 4;

        //������ �� ����

        app[0] = new mirror_reflection(Vect3f(_ks));
        app[1] = new diffuse_reflection(Vect3f(_kd));
        app[2] = new mirror_refraction(Vect3f(_kts), _n1, _n2);
        app[3] = new diffuse_refraction(Vect3f(_ktd), _n1, _n2);
        //printf("%f, %f, %f, %f\n", ks, kd, kts, ktd);

        //n = _n;
    }

    int diffusion(ray* ray, Vect3f normal) { //����������� - ���������� ��� �������
        double r = random.DRnd() * 1.2; //��������� � ����� ������
        int i;
        int j = -1;
        if (r > 1) {
            return -1;
        }
        for (i = 0; i < count_event; i++) {
            r = app[i] -> choose_event(r, ray);
            if (r <= 0) {
                j = app[i]->diffusion(ray, normal);
                break;
            }
        }
        //���������� (������ ����)
        return j; //�����, ������������ ��� ��������� ������, ������ ��� ������ �����������
    }

    Vect3f luminance(Vect3f E) { // ������� ������� +�������, ����������� ���������, ����������� ����������
        Vect3f res = Vect3f(rgb.x * E.x, rgb.y * E.y, rgb.z * E.z);
        res /= PI;
        return res;
    }
};

class camera { // ����� ������
public:
    Point3f position; // ����
    Point3f topleft; // ����� ������� ����
    int resolutionX, resolutionY; // ���������� �� � � �� �
    Vect3f dx, dy; // ������� �������

public:
    camera() {
        position = Point3f(5, 0, 0);
        topleft = Point3f(-5, -10, 10);

        resolutionX = 256;
        resolutionY = 256;

        dx = Vect3f(2 * -topleft.x / 256, 0, 0);
        dy = Vect3f(0, -2 * topleft.y / 256, 0);
    }
    camera(int resX, int resY, float anglez) {
        position = Point3f(10, 0, 0);
        topleft = Point3f(-5, -10, 10);

        resolutionX = resX;
        resolutionY = resY;

        float i = -2 * topleft.y / resolutionX;
        float j = -2 * topleft.z / resolutionY;

        dx = Vect3f(0, i, 0);
        dy = Vect3f(0, 0, j);

        Matrix43f tm(1); // ������� ������ ������ ����
        //tm.RotationX(anglex * PI / 180);
        /*tm.RotationY(angley * PI / 180);*/
        tm.RotationZ(anglez * PI / 180);
        tm.PointTransform(topleft);
        tm.PointTransform(position);
        tm.VectorTransform(dx);
        tm.VectorTransform(dy);
    }

    ray getRay(int x, int y) { // ��� ���������� ����� ���� � xy-������� ������
        int val = 1000;
        Vect3f randx = dx * random.DRnd();
        Vect3f randy = dy * random.DRnd();
        Point3f screen_point = topleft + x * dx + randx + y * dy + randy; // ������� �� ������
        //Point3f screen_point = topleft + x * dx  + y * dy;
        return ray(position, screen_point);
    }
};

struct point_and_normal { // ��������� � ������ �����������, �������� � ����������� � id �������
    Point3f point;
    Vect3f N;
    int id;
};

bool accur_equal(Point3f a, Point3f b, int accuracy) { // ��������� ����� � ��������� �� accuracy
    float x1 = round(a.x * accuracy) / accuracy;
    float y1 = round(a.y * accuracy) / accuracy;
    float z1 = round(a.z * accuracy) / accuracy;

    float x2 = round(b.x * accuracy) / accuracy;
    float y2 = round(b.y * accuracy) / accuracy;
    float z2 = round(b.z * accuracy) / accuracy;

    if (x1 == x2 && y1 == y2 && z1 == z2)
        return 1;

    return 0;
}

class geometry { // ����� ���������
public:
    RTCScene sc; // embree �����

public:
    geometry() {
        RTCDevice device = rtcNewDevice(NULL);
        if (device == NULL)
        {
            printf("Device error %d: cannot create device\n", rtcGetDeviceError(NULL));
            return;
        }
        rtcSetDeviceErrorFunction(device, DeviceErrorFunction, NULL);

        sc = rtcNewScene(device);
    }
    geometry(light* lights, int nl) { // �������� ���������
        RTCDevice device = rtcNewDevice(NULL);
        if (device == NULL)
        {
            printf("Device error %d: cannot create device\n", rtcGetDeviceError(NULL));
            return;
        }
        rtcSetDeviceErrorFunction(device, DeviceErrorFunction, NULL);

        sc = rtcNewScene(device);

        /*RTCGeometry sphere1 = CreateSphere(device, Point3f(1, -5, 10), 5, 5, Matrix43f(1, 1, 1));
        rtcAttachGeometry(sc, sphere1);

        RTCGeometry sphere2 = CreateSphere(device, Point3f(1, -2, 3), 2, 5, Matrix43f(1, 1, 1));
        rtcAttachGeometry(sc, sphere2);*/


        /*RTCGeometry obj = ReadOBJ(device, "C:/Color_Pyramid.obj");
        rtcAttachGeometry(sc, obj);*/

        for (int i = 0; i < nl; i++) { // ������������� ����� � �����
            RTCGeometry light = lights[i].lightGeom(device);
            rtcAttachGeometry(sc, light);
        }

        /*RTCGeometry obj = ReadOBJ(device, "D:\plant.obj");*/
        RTCGeometry obj = ReadOBJ(device, "C:/capybara.obj", 4);
        rtcAttachGeometry(sc, obj);

        printf("ok");

        Matrix43f tr = Matrix43f(1, 1, 1);
        tr.RotationX(0 * PI / 180);

        RTCGeometry box = CreateBox(device, Point3f(-30, -30, -12), Vect3f(60, 60, 10), tr);
        rtcAttachGeometry(sc, box);

        /*RTCGeometry sphere1 = CreateSphere(device, Point3f(2, -4, 1), 0.7, 5, Matrix43f(1, 1, 1));
        rtcAttachGeometry(sc, sphere1);*/

        RTCGeometry sphere1 = CreateSphere(device, Point3f(3, -4, 0), 0.5, 5, Matrix43f(1, 1, 1));
        rtcAttachGeometry(sc, sphere1);

        RTCGeometry sphere2 = CreateSphere(device, Point3f(3, -4, -2), 0.5, 5, Matrix43f(1, 1, 1));        
        rtcAttachGeometry(sc, sphere2);

        RTCGeometry sphere3 = CreateSphere(device, Point3f(3, 0, -2), 0.5, 5, Matrix43f(1, 1, 1));        
        rtcAttachGeometry(sc, sphere3);

        RTCGeometry sphere4 = CreateSphere(device, Point3f(3, -2, -2), 0.5, 5, Matrix43f(1, 1, 1));        
        rtcAttachGeometry(sc, sphere4);

        RTCGeometry sphere5 = CreateSphere(device, Point3f(-1, -3, -2), 0.5, 5, Matrix43f(1, 1, 1));
        rtcAttachGeometry(sc, sphere5);

        /*Matrix43f tr = Matrix43f(1, 1, 1);
        tr.RotationX(0 * PI / 180);

        RTCGeometry box = CreateBox(device, Point3f(5, 0, -1), Vect3f(0.3, 2, 2), tr);
        rtcAttachGeometry(sc, box);*/

        //RTCGeometry sphere3 = CreateSphere(device, Point3f(3, 0, 0), 1, 5, Matrix43f(1, 1, 1));
        //rtcAttachGeometry(sc, sphere3);        

        //RTCGeometry sphere4 = CreateSphere(device, Point3f(2, 3, -1), 1, 5, Matrix43f(1, 1, 1));
        //rtcAttachGeometry(sc, sphere4);

        rtcCommitScene(sc);
    }

    point_and_normal Intersection(ray r) { // ����� ����������� ���� � ����������
        struct RTCRayHit rayhit;
        rayhit.ray.org_x = r.start.x;
        rayhit.ray.org_y = r.start.y;
        rayhit.ray.org_z = r.start.z;
        rayhit.ray.dir_x = r.dir.x;
        rayhit.ray.dir_y = r.dir.y;
        rayhit.ray.dir_z = r.dir.z;
        rayhit.ray.tnear = 0;
        rayhit.ray.tfar = (float)MathF::MAX_VALUE;
        rayhit.ray.mask = 0;
        rayhit.ray.flags = 0;
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

        struct RTCIntersectContext context;
        rtcInitIntersectContext(&context);

        rtcIntersect1(sc, &context, &rayhit);
        Point3f point;
        Vect3f N;
        int id = -1;

        if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
            point = r.start + r.dir * (rayhit.ray.tfar);
            N = Vect3f(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z);

            if (rayhit.hit.geomID == 3) {
                N = (1 - rayhit.hit.u - rayhit.hit.v) * VN[normals[rayhit.hit.primID].x] +
                rayhit.hit.u * VN[normals[rayhit.hit.primID].y] +
                rayhit.hit.v * VN[normals[rayhit.hit.primID].z];
            }           

            N = N.Normalize();
            id = rayhit.hit.geomID;
            //rayhit.hit.
        }
        else {
            point = Point3f(0);
            N = Vect3f(0);
            id = -1;
        }

        return point_and_normal{ point, N, id };
    }

    bool is_not_shadow(Point3f p, ray light) { // ��� �� ���� �� ����� ��������� ����� � ���� �����
        struct RTCRayHit rayhit;
        rayhit.ray.org_x = light.start.x;
        rayhit.ray.org_y = light.start.y;
        rayhit.ray.org_z = light.start.z;
        rayhit.ray.dir_x = light.dir.x;
        rayhit.ray.dir_y = light.dir.y;
        rayhit.ray.dir_z = light.dir.z;
        rayhit.ray.tnear = 0;
        rayhit.ray.tfar = (float)MathF::MAX_VALUE;
        rayhit.ray.mask = 0;
        rayhit.ray.flags = 0;
        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

        struct RTCIntersectContext context;
        rtcInitIntersectContext(&context);

        rtcIntersect1(sc, &context, &rayhit);

        Point3f result_point = light.start + light.dir * (rayhit.ray.tfar);

        if (accur_equal(p, result_point, 1))
            return 1;
        else
            return 0;
    }
};

class pair {
public:
    int a, b;

public:
    pair() {
        a = 0;
        b = 0;
    }
    pair(int _a, int _b) {
        a = _a;
        b = _b;
    }
};

class scene { // �����-������� �����
public:
    Rnd random = Rnd();

    geometry geom;
    int light_n = 3;
    light* lights = new light[light_n];

    int op_n = 10;
    optprop* op = new optprop[op_n];
    camera cam;
    int rx = 1000, ry = 1000;
    double** pdf;
    double** prob;

    double* lightpdf;
    double* lightprob;

public:
    scene() {
        /*lights[0] = light(-15, -5, -7, 12800, 30000, 10000);
        lights[1] = light(20, 0, -7, 10000, 12800, 30000);
        lights[2] = light(20, -20, 30, 10000, 10000, 10000);*/
        /*lights[0] = light(3, 0, -7, 30000, 30000, 30000);
        lights[1] = light(-10, 0, -7, 10000, 10000, 10000);*/
        //lights[1] = light(1, -10, 10, 1840, 810, 2220);
        /*lights[0] = light(1, 10, 10, 10000, 10000, 10000);
        lights[1] = light(1, -10, 10, 10000, 10000, 10000);*/
        //lights[2] = light(-30, 0, 0, 10000, 10000, 10000);

       /* lights[0] = light(Point3f(20, -3, 10), Point3f(20, -5, 13), Point3f(20, -3, 11), 5000, 5000, 5000);
        lights[1] = light(Point3f(-2, -7, 11), Point3f(-2, -9, 12), Point3f(-4, -7, 13), 1000, 1000, 1000);
        lights[2] = light(Point3f(-7, -5, -6), Point3f(-6, -6, -6), Point3f(-5, -7, -7), 200, 200, 500);*/
        lights[0] = light(Point3f(4, -3, 1), Point3f(5, -4, 2), Point3f(6, -3, 1), 50000, 50000, 50000);        
        lights[1] = light(Point3f(-2, 0, 5), Point3f(-2, -2, 5), Point3f(-4, 0, 5), 5000, 10000, 5000);
        lights[2] = light(Point3f(-7, -5, 5), Point3f(-6, -6, 5), Point3f(-5, -7, 4), 50000, 10000, 10000);
        geom = geometry(lights, light_n);
        op[0] = optprop(5000, 1000, 1000, 1, 0, 0, 0, 0, 1, 1.5);
        op[1] = optprop(1000, 1000, 5000, 1, 0, 0, 0, 0, 1, 1.5);
        op[2] = optprop(1000, 1000, 5000, 1, 0, 0, 0, 0, 1, 1.5);

        op[3] = optprop(122, 62, 22, 0, 0, 0.1, 0, 0, 1, 1.5);
        op[4] = optprop(1, 1, 1.5, 0, 0, 0, 0.7, 0.1, 1, 1.9);
        op[5] = optprop(1, 1, 1, 0, 0, 0, 1, 0, 1, 1.5);

        op[6] = optprop(255, 79, 0, 0, 0.3, 0.3, 0, 0, 1, 1.5);
        op[7] = optprop(255, 79, 0, 0, 0.3, 0.3, 0, 0, 1, 1.5);
        op[8] = optprop(255, 79, 0, 0, 0.3, 0.3, 0, 0, 1, 1.5);
        op[9] = optprop(255, 79, 0, 0, 0.3, 0.3, 0, 0, 1, 1.5);
        /*op[5] = optprop(100, 255, 100, 0, 0.3, 0.1, 0, 0, 1, 1.5);
        op[6] = optprop(100, 100, 255, 0, 0.3, 0.1, 0, 0, 1, 1.5);*/
        /*op[1] = optprop(255, 100, 100, 0, 0);*/
        /*op[1] = optprop(100, 100, 255, 0);
        op[2] = optprop(255, 30, 255, 0);
        op[3] = optprop(100, 255, 100, 0);*/
        //op[0] = optprop(27, 50, 141);
        //op[1] = optprop(172, 143, 93);

        cam = camera(rx, ry, -45);

        pdf = new double* [rx];
        prob = new double* [rx];

        for (int i = 0; i < rx; i++) {
            pdf[i] = new double[ry];
            prob[i] = new double[ry + 1];

            pdf[i][0] = 1.0;
            prob[i][0] = 1.0;

            for (int j = 1; j < ry; j++) {
                pdf[i][j] = 1.0;
                prob[i][j] = prob[i][j - 1] + pdf[i][j];
            }

            prob[i][ry] = i == 0 ? prob[i][ry - 1] : prob[i][ry - 1] + prob[i - 1][ry];
        }

        for (int i = 0; i < rx; i++) {
            for (int j = 1; j < ry; j++) {
                prob[i][j] /= prob[i][ry - 1];
            }
            prob[i][ry] /= prob[rx - 1][ry];
        }

        lightpdf = new double[light_n];
        lightprob = new double[light_n];

        for (int i = 0; i < light_n; i++) {
            lightpdf[i] = lights[i].F();
            lightprob[i] = i == 0 ? lightpdf[i] : lightpdf[i] + lightprob[i - 1];
        }

        for (int i = 0; i < light_n; i++) {
            lightprob[i] /= lightprob[light_n - 1];
        }
        printf("hi\n");
    }

    void recount() {
        for (int i = 0; i < rx; i++) {
            prob[i][0] = pdf[i][0];

            for (int j = 1; j < ry; j++) {
                prob[i][j] = prob[i][j - 1] + pdf[i][j];
            }

            prob[i][ry] = i == 0 ? prob[i][ry - 1] : prob[i][ry - 1] + prob[i - 1][ry];
        }

        for (int i = 0; i < rx; i++) {
            for (int j = 1; j < ry; j++) {
                prob[i][j] /= prob[i][ry - 1];
            }
            prob[i][ry] /= prob[rx - 1][ry];
        }
    }

    int returnlight(double ir, int min, int max) {
        //printf("!a\n");
        int m = (max + min) / 2;

        if (m == 0 && ir <= lightprob[m])
            return m;
        else if (m == 0)
            return returnlight(ir, m + 1, max);

        if (ir > lightprob[m - 1] && ir <= lightprob[m])
            return m;

        if (ir > lightprob[m])
            return returnlight(ir, m + 1, max);
        else
            return returnlight(ir, min, m);
    }

    int returna(double ir, int min, int max) {
        //printf("!a\n");
        int m = (max + min) / 2;

        if (m == 0 && ir <= prob[m][ry])
            return m;
        else if (m == 0)
            return returna(ir, m + 1, max);

        if (ir > prob[m - 1][ry] && ir <= prob[m][ry])
            return m;

        if (ir > prob[m][ry])
            return returna(ir, m + 1, max);
        else
            return returna(ir, min, m);
    }

    int returnb(double jr, int a, int min, int max) {
        //printf("!b\n");
        int m = (max + min) / 2;

        if (m == 0 && jr <= prob[a][m])
            return m;
        else if (m == 0)
            return returnb(jr, a, m + 1, max);

        if (jr > prob[a][m - 1] && jr <= prob[a][m])
            return m;

        if (jr > prob[a][m])
            return returnb(jr, a, m + 1, max);
        else
            return returnb(jr, a, min, m);
    }

    pair getRandomRay() {
        double ir = random.DRnd();
        double jr = random.DRnd();

        //printf("%f, %f\n", ir, jr);

        int a = returna(ir, 0, rx - 1);
        int b = returnb(jr, a, 0, ry - 1);

        //printf("%d, %d\n", a, b);

        return pair(a, b);
    }

    int getRandomLight() {
        double ir = random.DRnd();

        int a = returnlight(ir, 0, light_n - 1);

        //printf("%d\n", a);

        return a;
    }
};

struct L_and_count {
    Vect3f L;
    int n;
};

// ����� ������ ����� ��� ������, ��������� �������

class render { // ����� �������
private:
    int NumResolution = 100;
    int NumLight = 1;

    TMatrix<Vect3f> m0; // ������ �������
    TMatrix<Vect3f> m1; // ������ ����
    TMatrix<Vect3f> m2; // ��������
    TMatrix<Vect3f> m3; // ��������� ���������

    TMatrix<Vect3f> sq0; // ������ �������
    TMatrix<Vect3f> sq1; // ������ ����
    TMatrix<Vect3f> sq2; // ��������
    TMatrix<Vect3f> sq3; // ��������� ���������

    TMatrix<int> num0;
    TMatrix<int> num1;
    TMatrix<int> num2;
    TMatrix<int> num3;

    TMatrix<int> ind;
    TMatrix<Vect3f> norm;

public:
    render() {
    }

    TMatrix<Vect3f> gauss_filter(TMatrix<Vect3f> m, int radius, int resx, int resy) {
        TMatrix<Vect3f> res = TMatrix<Vect3f>(resx, resy);

        for (int i = 0; i < resx; i++) {
            for (int j = 0; j < resy; j++) {
                if (ind[i][j] != -1) {
                    Vect3f sum = 0;
                    double w = 0;

                    Vect3f Lm = Vect3f(0);
                    int k = 0;

                    for (int in = i - radius; in < i + radius; in++) {
                        for (int jn = j - radius; jn < j + radius; jn++) {
                            if (in >= 0 && jn >= 0 && in < resx && jn < resy) {
                                k++;
                                Lm += m[in][jn];
                            }
                        }
                    }

                    Lm = Lm.Length() == 0 ? Vect3f(1) : Lm / k;

                    for (int in = i - radius; in < i + radius; in++) {
                        for (int jn = j - radius; jn < j + radius; jn++) {
                            if (in >= 0 && jn >= 0 && in < resx && jn < resy) {
                                if (ind[in][jn] == ind[i][j] && scalarMult(norm[in][jn], norm[i][j]) > 0.5) {
                                    double _w = Exp((-(in - i) * (in - i) - (jn - j) * (jn - j)) / 2 / radius / radius)
                                        * Exp(-radius) /** Exp((m[in][jn] - Lm).Length()/Lm.Length())*/;
                                    sum += m[in][jn] * _w;
                                    w += _w;
                                }
                            }
                        }
                    }
                    //res[i][j] = sum / 2 / PI / radius / radius;
                    res[i][j] = w == 0 ? sum : sum / w;
                }
                else {
                    res[i][j] = Vect3f(0);
                }
            }
        }
        return res;
    }

    int sigma(TMatrix<Vect3f> sq, TMatrix<Vect3f> m, TMatrix<int> num, int resx, int resy) {
        double res = 0;
        int k = 0;

        for (int i = 0; i < resx; i++) {
            for (int j = 0; j < resy; j++) {
                if (num[i][j] != 0) {
                    double S = (sq[i][j] / num[i][j] - (m[i][j] / num[i][j]) * (m[i][j] / num[i][j])).Length();
                    if (S < 0) S = 0;
                    res += Sqrt(S);
                    k++;
                }
            }
        }
        if (k == 0)
            k = 1;
        return Round(res / k);
    }

    Vect3f _L(scene sc, ray r, point_and_normal in, int a, int b) {
        geometry g = sc.geom;
        light* l = sc.lights;
        optprop* o = sc.op;

        int c = o[in.id].diffusion(&r, in.N);        

        if (c == -1) {
            //printf("hi0\n");
            return Vect3f(0);            
        }
        else {
            ray reflex = ray(in.point, r.dir);
            /*reflex.rgb = o[in.id].rgb;*/

            point_and_normal re = g.Intersection(reflex);

            if (re.id != -1) {
                //printf("popal\n");
                if (o[re.id].is_light) {
                    if (c == 0) {
                        return l[in.id].Lightness;
                    }
                    return Vect3f(0);
                }
                else {
                    Vect3f L = Vect3f(0);
                    int k = sc.getRandomLight();
                    for (int e = 0; e < NumLight; e++) {
                        Point3f rand = l[k].randomP();
                        ray li = ray(rand, re.point);

                        if (scalarMult(re.N, reflex.dir) * scalarMult(re.N, li.dir) >= 0) {
                            if (g.is_not_shadow(re.point, li)) {
                                L = o[re.id].luminance(l[k].illumination(re.point, rand, re.N));
                                L = k == 0 ? L / sc.lightprob[k] : L / (sc.lightprob[k] - sc.lightprob[k - 1]);

                                if (c%2 == 0) {
                                    m2[a][b] += L;
                                    sq2[a][b] += L*L;
                                    num2[a][b]++;
                                }
                                if (c%2 == 1) {
                                    m3[a][b] += L;
                                    sq3[a][b] += L * L;
                                    num3[a][b]++;
                                }
                                if (c == 2) {
                                    //printf("%d, %f %f %f\n", re.id, L.x, L.y, L.z);
                                }
                            }
                        }
                    }
                    /*if (c < 2)
                        return reflex.mirror_color(L + _L(sc, reflex, re, a, b));*/
                    return L + _L(sc, reflex, re, a, b);
                }
            }
        }
        return Vect3f(0);
    }

    void rendering(scene sc) { // ������� ��������
        printf("hi\n");
        geometry g = sc.geom;
        light* l = sc.lights;
        optprop* o = sc.op;
        camera c = sc.cam;

        TMatrix<Vect3f> m = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ��������

        m0 = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ������ �������
        m1 = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ������ ����
        m2 = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ��������
        m3 = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ��������� ���������

        sq0 = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ������ �������
        sq1 = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ������ ����
        sq2 = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ��������
        sq3 = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ��������� ���������

        ind = TMatrix<int>(c.resolutionX, c.resolutionY); // ������� ��������
        norm = TMatrix<Vect3f>(c.resolutionX, c.resolutionY); // ������� �����

        TMatrix<int> num = TMatrix<int>(c.resolutionX, c.resolutionY);

        num0 = TMatrix<int>(c.resolutionX, c.resolutionY);
        num1 = TMatrix<int>(c.resolutionX, c.resolutionY);
        num2 = TMatrix<int>(c.resolutionX, c.resolutionY);
        num3 = TMatrix<int>(c.resolutionX, c.resolutionY);

        for (int i = 0; i < c.resolutionX; i++) {
            for (int j = 0; j < c.resolutionY; j++) {
                m[i][j] = Vect3f(0);
                m0[i][j] = Vect3f(0);
                m1[i][j] = Vect3f(0);
                m2[i][j] = Vect3f(0);
                m3[i][j] = Vect3f(0);

                sq0[i][j] = Vect3f(0);
                sq1[i][j] = Vect3f(0);
                sq2[i][j] = Vect3f(0);
                sq3[i][j] = Vect3f(0);

                num[i][j] = 0;
                num0[i][j] = 0;
                num1[i][j] = 0;
                num2[i][j] = 0;
                num3[i][j] = 0;

                ind[i][j] = -1;
                norm[i][j] = Vect3f(0);
            }
        }

        int Num = c.resolutionX * c.resolutionY * NumResolution;

        for (int i = 0; i < Num; i++) {
            pair pair = sc.getRandomRay();

            ray r = c.getRay(pair.b, pair.a);

            /*printf("%d, %d\n", pair.a, pair.b);*/

            point_and_normal in = g.Intersection(r);

            if (in.id != -1) {
                ind[pair.a][pair.b] = in.id;
                norm[pair.a][pair.b] = in.N;
                //printf("popal\n");
                if (o[in.id].is_light) {
                    m[pair.a][pair.b] += l[in.id].Lightness;
                    num[pair.a][pair.b]++;
                    //printf("light\n");
                    m0[pair.a][pair.b] += l[in.id].Lightness;
                    sq0[pair.a][pair.b] += l[in.id].Lightness * l[in.id].Lightness;
                    num0[pair.a][pair.b]++;
                }
                else {
                    //printf("popal\n");
                    Vect3f L = Vect3f(0);
                    int k = sc.getRandomLight();
                    for (int e = 0; e < NumLight; e++) {
                        //printf("%d\n", e);
                        Point3f rand = l[k].randomP();
                        ray li = ray(rand, in.point);

                        if (scalarMult(in.N, r.dir) * scalarMult(in.N, li.dir) >= 0) {
                            if (g.is_not_shadow(in.point, li)) {
                                L = o[in.id].luminance(l[k].illumination(in.point, rand, in.N));
                                m1[pair.a][pair.b] += L;
                                sq1[pair.a][pair.b] += L*L;
                                num1[pair.a][pair.b]++;

                                L = k == 0 ? L / sc.lightprob[k] : L / (sc.lightprob[k] - sc.lightprob[k - 1]);
                                num[pair.a][pair.b]++;
                            }
                        }                      
                    }
                    L += _L(sc, r, in, pair.a, pair.b);
                    m[pair.a][pair.b] += L;
                    //num[pair.a][pair.b] += Lac.n;
                }
            }
            /*else {
                sc.pdf[pair.a][pair.b] = 0;
                sc.recount();
            }*/
        }

        for (int i = 0; i < c.resolutionX; i++) {
            for (int j = 0; j < c.resolutionY; j++) {
                m[i][j] = num[i][j] == 0 ? m[i][j] : m[i][j] / num[i][j];
                /*m0[i][j] = num0[i][j] == 0 ? m0[i][j] : m0[i][j] / num0[i][j];

                m1[i][j] = num1[i][j] == 0 ? m1[i][j] : m1[i][j] / num1[i][j];

                m2[i][j] = num2[i][j] == 0 ? m2[i][j] : m2[i][j] / num2[i][j];

                m3[i][j] = num3[i][j] == 0 ? m3[i][j] : m3[i][j] / num3[i][j];*/

            }
        }

        PathStr ps = PathStr("nit.nit"); // �������� �����
        WriteNITFile(ps, m, 0);
        PathStr ps0 = PathStr("nit0.nit"); // �������� �����
        WriteNITFile(ps0, m0, 0);

        m0 = gauss_filter(m0, 6/*sigma(sq0, m0, num0, c.resolutionX, c.resolutionY)*/, c.resolutionX, c.resolutionY);
        PathStr ps0_gauss = PathStr("nit0_gauss.nit"); // �������� �����
        WriteNITFile(ps0_gauss, m0, 0);

        PathStr ps1 = PathStr("nit1.nit"); // �������� �����
        WriteNITFile(ps1, m1, 0);

        m1 = gauss_filter(m1, 3/*sigma(sq1, m1, num1, c.resolutionX, c.resolutionY)*/, c.resolutionX, c.resolutionY);
        PathStr ps1_gauss = PathStr("nit1_gauss.nit"); // �������� �����
        WriteNITFile(ps1_gauss, m1, 0);

        PathStr ps2 = PathStr("nit2.nit"); // �������� �����
        WriteNITFile(ps2, m2, 0);
        
        m2 = gauss_filter(m2, 4/*sigma(sq2, m2, num2, c.resolutionX, c.resolutionY)*/, c.resolutionX, c.resolutionY);
        PathStr ps2_gauss = PathStr("nit2_gauss.nit"); // �������� �����
        WriteNITFile(ps2_gauss, m2, 0);

        PathStr ps3 = PathStr("nit3.nit"); // �������� �����
        WriteNITFile(ps3, m3, 0);

        m3 = gauss_filter(m3, 3, c.resolutionX, c.resolutionY);
        PathStr ps3_gauss = PathStr("nit3_gauss.nit"); // �������� �����
        WriteNITFile(ps3_gauss, m3, 0);

        for (int i = 0; i < c.resolutionX; i++) {
            for (int j = 0; j < c.resolutionY; j++) {
                if (num0[i][j] + num1[i][j] + num2[i][j] + num3[i][j] != 0)
                    m[i][j] = (m1[i][j] * num1[i][j] + m2[i][j] * num2[i][j] + m3[i][j] * num3[i][j] + m0[i][j] * num0[i][j]) / (num0[i][j] + num1[i][j] + num2[i][j] + num3[i][j]);
            }
        }

        PathStr ps_gauss = PathStr("nit_gauss.nit"); // �������� �����
        WriteNITFile(ps_gauss, m, 0);
    }
};

//////////////////////////////////////////////////////////////////////////
/// Program entry point.
int main()
{
    mem_init(NULL, NULL, "temp.mem");
    ev_init();
    col_init();

    render r = render();
    scene s = scene();
    r.rendering(s);

    col_term();
    ev_term();
    mem_close();

    return 0;
}
