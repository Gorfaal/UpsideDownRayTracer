#include <iostream>
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <stdlib.h>
#include <math.h>
#include <vtkRectilinearGrid.h>

using namespace std;

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
}

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
}

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
}

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
}

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);
}

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
     idx[0] = cellId%(dims[0]-1);
     idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
     idx[2] = cellId/((dims[0]-1)*(dims[1]-1));
}

struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};

struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins

    // Take in a value and applies the transfer function.
    // Step #1: figure out which bin "value" lies in.
    // If "min" is 2 and "max" is 4, and there are 10 bins, then
    //   bin 0 = 2->2.2
    //   bin 1 = 2.2->2.4
    //   bin 2 = 2.4->2.6
    //   bin 3 = 2.6->2.8
    //   bin 4 = 2.8->3.0
    //   bin 5 = 3.0->3.2
    //   bin 6 = 3.2->3.4
    //   bin 7 = 3.4->3.6
    //   bin 8 = 3.6->3.8
    //   bin 9 = 3.8->4.0
    // and, for example, a "value" of 3.15 would return the color in bin 5
    // and the opacity at "opacities[5]".
    int GetBin(double value)
    {
        double difference = max - min;
        double step = difference / numBins;

        for(int i = 0; i < numBins; i++)
        {
            if(value < (i * step) + min)
            {
                if(i != 0) return i - 1; else return 0;
            }
        }
    }
    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        int bin = GetBin(value);
        RGB[0] = colors[3*bin+0];
        RGB[1] = colors[3*bin+1];
        RGB[2] = colors[3*bin+2];
        opacity = opacities[bin];
    }
};

TransferFunction SetupTransferFunction(void)
{
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = {
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0,
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76
       };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins-1;
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0 ; k < 3 ; k++)
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }

    return rv;
}

Camera SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = 1;
    rv.up[2] = 0;
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;

    return rv;
}

void WriteImage(vtkImageData *img, const char *filename)
{
   string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInput(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

vtkImageData* NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    image->SetNumberOfScalarComponents(3);
    image->SetScalarType(VTK_UNSIGNED_CHAR);
    image->AllocateScalars();
    return image;
}

double* crossProduct(double* a, double *b)
{
    double *result = new double[3];
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
    return result;
}

double* vectorSubtract(double* a, double *b)
{
    double *result = new double[3];
    for(int i = 0; i < 3; i++)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

double* vectorTimes(double* a, double *b)
{
    double *result = new double[3];
    for(int i = 0; i < 3; i++)
    {
        result[i] = a[i] * b[i];
    }
    return result;
}

double vectorMagnitude(double* a)
{
    double sumOfSquares = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    double result = sqrt(sumOfSquares);
    return result;
}

double* vectorNormalize(double* a, double norm)
{
    double *result = new double[3];
    for(int i = 0; i < 3; i++)
    {
        result[i] = a[i] / norm;
    }
    return result;
}

double* vectorTimesDouble(double* a, double b)
{
    double *result = new double[3];
    for(int i = 0; i < 3; i++)
    {
        result[i] = a[i] * b;
    }
    return result;
}

double* vectorAddDouble(double* a, double b)
{
    double *result = new double[3];
    for(int i = 0; i < 3; i++)
    {
        result[i] = a[i] + b;
    }
    return result;
}

double* vectorAdd(double* a, double *b)
{
    double *result = new double[3];
    for(int i = 0; i < 3; i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}

float interpolate(float x, float v1, float v2){
    return v1 * (1 - x) + v2 * x;
}

double EvaluateFieldAtLocation(const double *pt, const int *dims, const float *X, const float *Y, const float *Z, const float *F)
{
    int xMin = -1;
    int yMin = -1;
    int zMin = -1;

    for(int i = 0; i < dims[0]; i++)
    {
        if(pt[0] <= X[i])
        {
            xMin = i - 1;
            break;
        }
    }

    for(int j = 0; j < dims[1]; j++)
    {
        if(pt[1] <= Y[j])
        {
            yMin = j - 1;
            break;
        }
    }

    for(int k = 0; k < dims[2]; k++)
    {
        if(pt[2] <= Z[k])
        {
            zMin = k- 1;
            break;
        }
    }

    if(xMin < 0 || yMin < 0 || zMin < 0)
    {
        return 0; //assume out of range so just return 0
    }

    //FOLLOWING WIKIPEDIA NOTATION/ALGORITHM
    //http://en.wikipedia.org/wiki/Trilinear_interpolation

    int c000LogicalIndex[3];
    c000LogicalIndex[0] = xMin;
    c000LogicalIndex[1] = yMin;
    c000LogicalIndex[2] = zMin;
    int c000Index = GetPointIndex(c000LogicalIndex, dims);
    double c000Value = F[c000Index];

    int c001LogicalIndex [3];
    c001LogicalIndex[0] = xMin;
    c001LogicalIndex[1] = yMin;
    c001LogicalIndex[2] = zMin + 1;
    int c001Index = GetPointIndex(c001LogicalIndex, dims);
    double c001Value = F[c001Index];

    int c010LogicalIndex [3];
    c010LogicalIndex[0] = xMin;
    c010LogicalIndex[1] = yMin + 1;
    c010LogicalIndex[2] = zMin;
    int c010Index = GetPointIndex(c010LogicalIndex, dims);
    double c010Value = F[c010Index];

    int c011LogicalIndex [3];
    c011LogicalIndex[0] = xMin;
    c011LogicalIndex[1] = yMin + 1;
    c011LogicalIndex[2] = zMin + 1;
    int c011Index = GetPointIndex(c011LogicalIndex, dims);
    double c011Value = F[c011Index];

    int c100LogicalIndex [3];
    c100LogicalIndex[0] = xMin + 1;
    c100LogicalIndex[1] = yMin;
    c100LogicalIndex[2] = zMin;
    int c100Index = GetPointIndex(c100LogicalIndex, dims);
    double c100Value = F[c100Index];

    int c101LogicalIndex [3];
    c101LogicalIndex[0] = xMin + 1;
    c101LogicalIndex[1] = yMin;
    c101LogicalIndex[2] = zMin + 1;
    int c101Index = GetPointIndex(c101LogicalIndex, dims);
    double c101Value = F[c101Index];

    int c110LogicalIndex [3];
    c110LogicalIndex[0] = xMin + 1;
    c110LogicalIndex[1] = yMin + 1;
    c110LogicalIndex[2] = zMin;
    int c110Index = GetPointIndex(c110LogicalIndex, dims);
    double c110Value = F[c110Index];

    int c111LogicalIndex [3];
    c111LogicalIndex[0] = xMin + 1;
    c111LogicalIndex[1] = yMin + 1;
    c111LogicalIndex[2] = zMin + 1;
    int c111Index = GetPointIndex(c111LogicalIndex, dims);
    double c111Value = F[c111Index];

    double x = pt[0];
    double y = pt[1];
    double z = pt[2];

    double fx = (x - X[xMin]) / (X[xMin + 1] - X[xMin]);
    double fy = (y - Y[yMin]) / (Y[yMin + 1] - Y[yMin]);
    double fz = (z - Z[zMin]) / (Z[zMin + 1] - Z[zMin]);

    double c00 = interpolate(fx, c000Value, c100Value);
    double c01 = interpolate(fx, c001Value, c101Value);
    double c10 = interpolate(fx, c010Value, c110Value);
    double c11 = interpolate(fx, c011Value, c111Value);

    double c0 = interpolate(fy, c00, c10);
    double c1 = interpolate(fy, c01, c11);

    double c = interpolate(fz, c0, c1);
    return c;
}

int main()
{
    TransferFunction tf = SetupTransferFunction();
    Camera camera = SetupCamera();
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("astro64.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    int height = 300;
    int width = 300;
    double sampleSize = 250;
    vtkImageData *image = NewImage(width, height);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = height * width;


    for (int i = 0 ; i < npixels*3; i++)
    {
        buffer[i] = 0;
    }


    for(int px = 0; px < width; px++)
    {
        for(int py = 0; py < height; py++)
        {
            int currentPixel = 3 * (py * width + px);
            double *look = vectorSubtract(camera.focus, camera.position);
            double *lookCrossUp = crossProduct(look, camera.up);
            double magLookCrossUp = vectorMagnitude(lookCrossUp);
            double *ru = vectorNormalize(lookCrossUp, magLookCrossUp);
            double *lookCrossRu = crossProduct(look, ru);
            double magLookCrossRu = vectorMagnitude(lookCrossRu);
            double *rv = vectorNormalize(lookCrossRu, magLookCrossRu);
            double magLook = vectorMagnitude(look);

            double *lookNorm = vectorNormalize(look, magLook);
            double angle_rad = camera.angle/360.0*(2*3.14159);
            double *rx = vectorTimesDouble(ru,(2 * tan(angle_rad / 2) / width));
            double *ry = vectorTimesDouble(rv,(2 * tan(angle_rad / 2) / height));

            double *pixelXCompenent= vectorTimesDouble(rx,(2 * px + 1 - width) / 2);
            double *pixelYCompenent= vectorTimesDouble(ry, (2 * py + 1 - height) / 2);

            double *a = vectorAdd(pixelYCompenent, pixelXCompenent);
            double *ray = vectorAdd(a, lookNorm);

            double distance = camera.far - camera.near;
            double sampleStep = distance / sampleSize;

            unsigned char rayRGB[3];
            rayRGB[0] = 0;
            rayRGB[1] = 0;
            rayRGB[2] = 0;
            double rayOpacity = 0;
            for(double s = camera.near; s <= camera.far; s+= sampleStep)
            {
                double *currentPosition = vectorTimesDouble(ray, s);
                double *updatedPosition = vectorAdd(camera.position, currentPosition);
                double val = EvaluateFieldAtLocation(updatedPosition, dims, X, Y, Z, F);

                unsigned char sampleRGB[3];
                double sampleOpacity = 0;
                tf.ApplyTransferFunction(val, sampleRGB, sampleOpacity);
                if(sampleOpacity != 0)
                {

                    unsigned char newRed = (rayOpacity*rayRGB[0] + (1-rayOpacity))*sampleRGB[0];
                    unsigned char newGreen = (rayOpacity*rayRGB[1] + (1-rayOpacity))*sampleRGB[1];
                    unsigned char newBlue = (rayOpacity*rayRGB[2] + (1-rayOpacity))*sampleRGB[2];
                    rayOpacity = rayOpacity + (1 - rayOpacity) * sampleOpacity;
                    rayRGB[0] = newRed;
                    rayRGB[1] = newGreen;
                    rayRGB[2] = newBlue;


                }
                delete [] currentPosition;
                delete [] updatedPosition;

                if(rayOpacity >= 1)
                {
                    break;
                }
            }

            buffer[currentPixel] = rayRGB[0];
            buffer[currentPixel + 1] = rayRGB[1];
            buffer[currentPixel + 2] = rayRGB[2];
            delete [] look;
            delete [] lookCrossUp;
            delete [] ru;
            delete [] lookCrossRu;
            delete [] rv;
            delete [] lookNorm;
            delete [] rx;
            delete [] ry;
            delete [] pixelXCompenent;
            delete [] pixelYCompenent;
            delete [] a;
            delete [] ray;
        }
    }

    WriteImage(image, "astro64");
}

