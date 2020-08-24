/*
 * Copyright (C) 2017-2020 HERE Europe B.V.
 * Licensed under Apache 2.0, see full license in LICENSE
 * SPDX-License-Identifier: Apache-2.0
 */
import {
    EarthConstants,
    GeoCoordinates,
    GeoPolygon,
    GeoPolygonCoordinates,
    Projection,
    ProjectionType
} from "@here/harp-geoutils";
import { assert, LoggerManager, Math2D } from "@here/harp-utils";
import {
    EllipseCurve,
    Frustum,
    Line3,
    MathUtils,
    Matrix4,
    PerspectiveCamera,
    Plane,
    Ray,
    Vector2,
    Vector3
} from "three";

import { TileCorners } from "./geometry/TileGeometryCreator";
import { MapViewUtils } from "./Utils";

function computeLongitudeSpan(geoStart: GeoCoordinates, geoEnd: GeoCoordinates): number {
    // TODO: Probably will need fix, checking camera orientation, frustum side and winding.

    const minLongitude = Math.min(geoStart.longitude, geoEnd.longitude);
    const maxLongitude = Math.max(geoStart.longitude, geoEnd.longitude);

    return Math.min(maxLongitude - minLongitude, 360 + minLongitude - maxLongitude);
}

function computeDivisionCount(geoStart: GeoCoordinates, geoEnd: GeoCoordinates): number {
    // TODO: Refine subdivision criteria.

    const maxLatitudeSpan = 20;
    const maxLongitudeSpan = 5;

    const latitudeSpan = Math.abs(geoEnd.latitude - geoStart.latitude);
    const longitudeSpan = computeLongitudeSpan(geoStart, geoEnd);
    return Math.ceil(Math.max(latitudeSpan / maxLatitudeSpan, longitudeSpan / maxLongitudeSpan));
}

const logger = LoggerManager.instance.create("BoundsGenerator");

// Keep ccw order!
enum Side {
    Bottom,
    Right,
    Top,
    Left
}

function nextSide(side: Side): Side {
    return (side + 1) % 4;
}

function previousSide(side: Side): Side {
    return (side + 3) % 4;
}

const ccwCanvasCornersNDC: Array<{ x: number; y: number }> = [
    { x: -1, y: -1 }, // bottom left
    { x: 1, y: -1 }, // bottom right
    { x: 1, y: 1 }, // top right
    { x: -1, y: 1 } // top left
];

class SphereHorizon {
    private readonly m_matrix: Matrix4;
    private readonly m_horizonCircle: EllipseCurve;
    private readonly m_normalToTangentAngle: number;
    private readonly m_distanceToHorizonCenter: number;
    private readonly m_intersections: number[][] = [];
    private m_cameraPitch?: number;
    private m_hFovVertical?: number;
    private m_hFovHorizontal?: number;

    constructor(
        private readonly m_camera: PerspectiveCamera,
        private readonly m_projection: Projection
    ) {
        const earthRadiusSq = EarthConstants.EQUATORIAL_RADIUS * EarthConstants.EQUATORIAL_RADIUS;
        const xAxis = new Vector3().setFromMatrixColumn(m_camera.matrixWorld, 0).normalize();
        const zAxis = m_camera.position.clone().normalize();
        const yAxis = new Vector3().crossVectors(zAxis, xAxis);

        const cameraHeight = m_camera.position.length();
        this.m_normalToTangentAngle = Math.asin(EarthConstants.EQUATORIAL_RADIUS / cameraHeight);

        const tangentDistance = Math.sqrt(cameraHeight * cameraHeight - earthRadiusSq);
        this.m_distanceToHorizonCenter = tangentDistance * Math.cos(this.m_normalToTangentAngle);
        const horizonCenterLength = cameraHeight - this.m_distanceToHorizonCenter;
        const horizonRadius = Math.sqrt(earthRadiusSq - horizonCenterLength * horizonCenterLength);
        const horizonCenter = new Vector3().copy(zAxis).setLength(horizonCenterLength);

        this.m_matrix = new Matrix4().makeBasis(xAxis, yAxis, zAxis).setPosition(horizonCenter);
        this.m_horizonCircle = new EllipseCurve(
            0,
            0,
            horizonRadius,
            horizonRadius,
            0,
            Math.PI * 2,
            false,
            0
        );

        this.computeIntersections();
    }

    getPoint(t: number, target: Vector3 = new Vector3()): Vector3 {
        const point2d = this.m_horizonCircle.getPoint(t, new Vector2());

        target.set(point2d.x, point2d.y, 0);
        target.applyMatrix4(this.m_matrix);
        return target;
    }

    getDivisionPoints(callback: (point: Vector3) => void, tStart = 0, tEnd = 1) {
        const divisionCount = Math.ceil((tEnd - tStart) / 0.05); // 20 for a whole horizon circle.
        this.m_horizonCircle.aStartAngle = tStart * Math.PI * 2;
        this.m_horizonCircle.aEndAngle = tEnd * Math.PI * 2;
        const points = this.m_horizonCircle.getPoints(divisionCount);
        --points.length; // discard end point.
        const point3d = new Vector3();
        points.forEach((point2d: Vector2) => {
            point3d.set(point2d.x, point2d.y, 0);
            point3d.applyMatrix4(this.m_matrix);
            callback(point3d);
        });

        this.m_horizonCircle.aStartAngle = 0;
        this.m_horizonCircle.aEndAngle = 2 * Math.PI;
    }

    isFullyVisible(): boolean {
        // Only need to check either left or right.
        return (
            this.isTangentVisible(Side.Top) &&
            this.isTangentVisible(Side.Bottom) &&
            this.isTangentVisible(Side.Left)
        );
    }

    getSideIntersections(side: Side): number[] {
        return this.m_intersections[side];
    }

    private isTangentVisible(side: Side): boolean {
        if (side === Side.Top || side === Side.Bottom) {
            const eyeToTangentAngle =
                side === Side.Top
                    ? this.m_normalToTangentAngle - this.cameraPitch
                    : this.m_normalToTangentAngle + this.cameraPitch;
            return this.hFovVertical >= Math.abs(eyeToTangentAngle);
        } else {
            const eyeToTangentAngle = this.m_normalToTangentAngle;
            return (
                this.hFovHorizontal >= Math.abs(eyeToTangentAngle) &&
                this.cameraPitch <= this.hFovVertical
            );
        }
    }

    private getTangentOnSide(side: Side): number {
        switch (side) {
            case Side.Bottom:
                return 0.75;
            case Side.Right:
                return 0;
            case Side.Top:
                return 0.25;
            case Side.Left:
                return 0.5;
        }
    }

    private computeIntersections() {
        const rSq = this.m_horizonCircle.xRadius * this.m_horizonCircle.xRadius;
        const yBottom =
            this.m_distanceToHorizonCenter * Math.tan(this.cameraPitch - this.hFovVertical);
        let tTR: number | undefined;
        let tBR: number | undefined;

        for (let side = Side.Bottom; side < 4; side++) {
            const sideIntersections: number[] = [];
            if (this.isTangentVisible(side)) {
                sideIntersections.push(this.getTangentOnSide(side));
            } else {
                switch (side) {
                    case Side.Bottom: {
                        const x = Math.sqrt(rSq - yBottom * yBottom);
                        const t = Math.atan2(yBottom, x) / (Math.PI * 2);
                        sideIntersections.push(0.5 - t, t);
                        break;
                    }
                    case Side.Right: {
                        const eyeToHorizon =
                            this.m_distanceToHorizonCenter / Math.cos(this.cameraPitch);
                        const yR = this.m_distanceToHorizonCenter * Math.tan(this.cameraPitch);
                        const xR = eyeToHorizon * Math.tan(this.hFovHorizontal);
                        const dB =
                            (this.m_distanceToHorizonCenter * Math.cos(this.hFovVertical)) /
                            Math.cos(this.cameraPitch - this.hFovVertical);
                        const xBr = (xR * dB) / eyeToHorizon;
                        const yBr = yBottom;
                        // line: (xR,yR) - (xBr, yBr).
                        const intersections = Math2D.intersectLineAndCircle(
                            xBr,
                            yBr,
                            xR,
                            yR,
                            this.m_horizonCircle.xRadius
                        );
                        const yTopRight = intersections!.y1;
                        if (-yTopRight >= yBottom) {
                            tBR =
                                Math.atan2(intersections!.y2!, intersections!.x2!) / (Math.PI * 2);
                            sideIntersections.push(tBR);
                        }
                        tTR = Math.atan2(intersections!.y1, intersections!.x1) / (Math.PI * 2);
                        sideIntersections.push(tTR);
                        break;
                    }
                    case Side.Top: {
                        const yTop =
                            this.m_distanceToHorizonCenter *
                            Math.tan(this.cameraPitch + this.hFovVertical);
                        const x = Math.sqrt(rSq - yTop * yTop);
                        const t = Math.atan2(yTop, x) / (Math.PI * 2);
                        sideIntersections.push(t, 0.5 - t);
                        break;
                    }
                    case Side.Left: {
                        sideIntersections.push(0.5 - tTR!);
                        if (tBR !== undefined) {
                            sideIntersections.push(0.5 - tBR);
                        }
                        break;
                    }
                }
            }

            this.m_intersections.push(sideIntersections);
        }
    }

    private get cameraPitch(): number {
        if (this.m_cameraPitch === undefined) {
            this.m_cameraPitch = MapViewUtils.extractAttitude(
                { projection: this.m_projection },
                this.m_camera
            ).pitch;
        }
        return this.m_cameraPitch;
    }

    private get hFovVertical(): number {
        if (this.m_hFovVertical === undefined) {
            this.m_hFovVertical = MathUtils.degToRad(this.m_camera.fov / 2);
        }
        return this.m_hFovVertical;
    }

    private get hFovHorizontal(): number {
        if (this.m_hFovHorizontal === undefined) {
            this.m_hFovHorizontal =
                MapViewUtils.calculateHorizontalFovByVerticalFov(
                    this.hFovVertical * 2,
                    this.m_camera.aspect
                ) / 2;
        }
        return this.m_hFovHorizontal;
    }
}

/**
 * Generates Bounds for a camera view and a projection
 *
 * @beta, @internal
 */
export class BoundsGenerator {
    private readonly m_groundPlaneNormal = new Vector3(0, 0, 1);
    private readonly m_groundPlane = new Plane(this.m_groundPlaneNormal.clone());

    constructor(
        private readonly m_camera: PerspectiveCamera,
        private m_projection: Projection,
        public tileWrappingEnabled: boolean = false
    ) {}

    set projection(projection: Projection) {
        this.m_projection = projection;
    }

    /**
     * Generates an Array of GeoCoordinates covering the visible map.
     * The coordinates are sorted to ccw winding, so a polygon could be drawn with them.
     */
    generate(): GeoPolygon | undefined {
        return this.m_projection.type === ProjectionType.Planar
            ? this.generateOnPlane()
            : this.generateOnSphere();
    }

    private createPolygon(
        coordinates: GeoCoordinates[],
        sort: boolean,
        normalize: boolean = false
    ): GeoPolygon | undefined {
        if (coordinates.length > 2) {
            if (normalize) {
                this.normalizeCoordinates(coordinates);
            }
            return new GeoPolygon(coordinates as GeoPolygonCoordinates, sort);
        }
        return undefined;
    }

    private normalizeCoordinates(coordinates: GeoCoordinates[]) {
        // let minLongitude = Infinity;
        // let maxLongitude = -Infinity;
        // for (const geoCoords of coordinates) {
        //     minLongitude = Math.min(minLongitude, geoCoords.longitude);
        //     maxLongitude = Math.max(maxLongitude, geoCoords.longitude);
        // }
        // if (minLongitude >= 0 || maxLongitude < 0) {
        //     return;
        // }
        // const cameraTarget = MapViewUtils.rayCastWorldCoordinates(
        //     { camera: this.m_camera, projection: this.m_projection },
        //     0,
        //     0
        // );
        // if (cameraTarget) {
        //     const targetLongitude = this.m_projection.unprojectPoint(cameraTarget).longitude;
        //     if (targetLongitude > -90 && targetLongitude < 90) {
        //         return;
        //     }
        // } else if (maxLongitude - minLongitude < 360 + minLongitude - maxLongitude) {
        //     return;
        // }
        // for (const geoCoord of coordinates) {
        //     if (geoCoord.longitude < 0) {
        //         geoCoord.longitude = 360 + geoCoord.longitude;
        //     }
        // }

        const antimerCrossIndex = coordinates.findIndex((val: GeoCoordinates, index: number) => {
            const prevLonIndex = index === 0 ? coordinates.length - 1 : index - 1;
            const prevLon = coordinates[prevLonIndex].longitude;
            const lon = val.longitude;

            if (prevLon > 90 && lon < -90) {
                return true;
            }
        });
        if (antimerCrossIndex < 0) {
            return;
        }

        for (let i = 0; i < coordinates.length; i++) {
            const index = (antimerCrossIndex + i) % coordinates.length;
            const currentLon = coordinates[index].longitude;
            coordinates[index].longitude += 360;
            const nextLon = coordinates[(index + 1) % coordinates.length].longitude;

            if (currentLon < -90 && nextLon > 90) {
                // new crossing in opposite direction, stop.
                break;
            }
        }
    }

    private addSideSegmentSubdivisions(
        coordinates: GeoCoordinates[],
        NDCStart: { x: number; y: number },
        NDCEnd: { x: number; y: number },
        geoStart: GeoCoordinates,
        geoEnd: GeoCoordinates
    ) {
        coordinates.push(geoStart);

        const divisionCount = computeDivisionCount(geoStart, geoEnd);
        if (divisionCount <= 1) {
            return;
        }

        const NDCStep = new Vector2(NDCEnd.x - NDCStart.x, NDCEnd.y - NDCStart.y).multiplyScalar(
            1 / divisionCount
        );

        const NDCDivision = new Vector2(NDCStart.x, NDCStart.y);
        for (let i = 0; i < divisionCount - 1; i++) {
            NDCDivision.add(NDCStep);
            const intersection = MapViewUtils.rayCastWorldCoordinates(
                { camera: this.m_camera, projection: this.m_projection },
                NDCDivision.x,
                NDCDivision.y
            );
            if (intersection) {
                coordinates.push(this.m_projection.unprojectPoint(intersection));
            }
        }
    }

    private addSideIntersectionsOnSphere(
        coordinates: GeoCoordinates[],
        side: Side,
        geoStartCorner?: GeoCoordinates,
        geoEndCorner?: GeoCoordinates,
        horizon?: SphereHorizon
    ) {
        const startNDCCorner = ccwCanvasCornersNDC[side];
        const endNDCCorner = ccwCanvasCornersNDC[nextSide(side)];

        if (geoStartCorner && geoEndCorner) {
            this.addSideSegmentSubdivisions(
                coordinates,
                startNDCCorner,
                endNDCCorner,
                geoStartCorner,
                geoEndCorner
            );
            return;
        }

        if (!horizon) {
            return;
        }

        const horizonIntersections = horizon.getSideIntersections(side);
        if (horizonIntersections.length === 0) {
            return;
        }

        if (geoStartCorner) {
            // There should only be one horizon intersection, if there's 2 take the last one.
            const worldHorizonPoint = horizon.getPoint(
                horizonIntersections[horizonIntersections.length - 1]
            );
            const geoHorizonPoint = this.m_projection.unprojectPoint(worldHorizonPoint);
            this.addSideSegmentSubdivisions(
                coordinates,
                startNDCCorner,
                worldHorizonPoint.project(this.m_camera),
                geoStartCorner,
                geoHorizonPoint
            );
        } else {
            // Subdivide horizon from previous side last horizon intersection to this side first's.
            const prevSideIntersections = horizon.getSideIntersections(previousSide(side));
            assert(prevSideIntersections.length > 0);

            horizon.getDivisionPoints(
                point => {
                    coordinates.push(this.m_projection.unprojectPoint(point));
                },
                prevSideIntersections[prevSideIntersections.length - 1],
                horizonIntersections[0]
            );
        }

        if (horizonIntersections.length > 1) {
            // Subdivide side segment between two horizon intersections.
            const worldHorizonStart = horizon.getPoint(horizonIntersections[0]);
            const worldHorizonEnd = horizon.getPoint(horizonIntersections[1]);
            const geoHorizonStart = this.m_projection.unprojectPoint(worldHorizonStart);
            const geoHorizonEnd = this.m_projection.unprojectPoint(worldHorizonEnd);

            this.addSideSegmentSubdivisions(
                coordinates,
                worldHorizonStart.project(this.m_camera),
                worldHorizonEnd.project(this.m_camera),
                geoHorizonStart,
                geoHorizonEnd
            );
        }

        if (geoEndCorner) {
            const worldHorizonPoint = horizon.getPoint(horizonIntersections[0]);
            const geoHorizonPoint = this.m_projection.unprojectPoint(worldHorizonPoint);
            this.addSideSegmentSubdivisions(
                coordinates,
                worldHorizonPoint.project(this.m_camera),
                endNDCCorner,
                geoHorizonPoint,
                geoEndCorner
            );
        }
    }

    private findBoundsIntersectionsOnSphere(): GeoCoordinates[] {
        assert(this.m_projection.type === ProjectionType.Spherical);

        const cornerCoordinates: GeoCoordinates[] = [];
        const coordinates: GeoCoordinates[] = [];

        this.addCanvasCornerIntersection(cornerCoordinates);

        // Asumptions: No camera roll. Either 0, 2 (bl,br) or 4 corners intersect.

        const horizon =
            cornerCoordinates.length < 4
                ? new SphereHorizon(this.m_camera, this.m_projection)
                : undefined;

        if (cornerCoordinates.length === 0 && horizon!.isFullyVisible()) {
            horizon!.getDivisionPoints(point => {
                coordinates.push(this.m_projection.unprojectPoint(point));
            });
            return coordinates;
        }

        cornerCoordinates.length = 4;
        for (let side = Side.Bottom; side < 4; side++) {
            const startCorner = cornerCoordinates[side];
            const endCorner = cornerCoordinates[nextSide(side)];
            this.addSideIntersectionsOnSphere(coordinates, side, startCorner, endCorner, horizon);
        }
        return coordinates;
    }

    private generateOnSphere(): GeoPolygon | undefined {
        assert(this.m_projection.type === ProjectionType.Spherical);

        const coordinates = this.findBoundsIntersectionsOnSphere();
        return this.createPolygon(coordinates, false, true);
    }

    private generateOnPlane(): GeoPolygon | undefined {
        //!!!!!!!ALTITUDE IS NOT TAKEN INTO ACCOUNT!!!!!!!!!
        const coordinates: GeoCoordinates[] = [];

        // 1.) Raycast into all four corners of the canvas
        //     => if an intersection is found, add it to the polygon
        this.addCanvasCornerIntersection(coordinates);

        // => All 4 corners found an intersection, therefore the screen is covered with the map
        // and the polygon complete
        if (coordinates.length === 4) {
            return this.createPolygon(coordinates, true);
        }

        //2.) Raycast into the two corners of the horizon cutting the canvas sides
        //    => if an intersection is found, add it to the polygon
        this.addHorizonIntersection(coordinates);

        //Setup the frustum for further checks
        const frustum = new Frustum().setFromProjectionMatrix(
            new Matrix4().multiplyMatrices(
                this.m_camera.projectionMatrix,
                this.m_camera.matrixWorldInverse
            )
        );

        // Setup the world corners for further checks.
        // Cast to TileCorners as it cannot be undefined here, due to the forced
        // PlanarProjection above
        const worldCorners: TileCorners = this.getWorldConers(this.m_projection) as TileCorners;

        if (!this.tileWrappingEnabled) {
            // 3.) If no wrapping, check if any corners of the world plane are inside the view
            //     => if true, add it to the polygon
            [worldCorners.ne, worldCorners.nw, worldCorners.se, worldCorners.sw].forEach(corner => {
                this.addPointInFrustum(corner, frustum, coordinates);
            });
        }

        //4.) Check for any edges of the world plane intersecting with the frustum?
        //    => if true, add to polygon

        if (!this.tileWrappingEnabled) {
            // if no tile wrapping:
            //       check with limited lines around the world edges
            [
                new Line3(worldCorners.sw, worldCorners.se), // south edge
                new Line3(worldCorners.ne, worldCorners.nw), // north edge
                new Line3(worldCorners.se, worldCorners.ne), // east edge
                new Line3(worldCorners.nw, worldCorners.sw) //  west edge
            ].forEach(edge => {
                this.addFrustumIntersection(edge, frustum, coordinates);
            });
        } else {
            // if tile wrapping:
            //       check for intersections with rays along the south and north edges
            const directionEast = new Vector3() //west -> east
                .subVectors(worldCorners.sw, worldCorners.se)
                .normalize();
            const directionWest = new Vector3() //east -> west
                .subVectors(worldCorners.se, worldCorners.sw)
                .normalize();

            [
                new Ray(worldCorners.se, directionEast), // south east ray
                new Ray(worldCorners.se, directionWest), // south west ray
                new Ray(worldCorners.ne, directionEast), // north east ray
                new Ray(worldCorners.ne, directionWest) //  north west ray
            ].forEach(ray => {
                this.addFrustumIntersection(ray, frustum, coordinates);
            });
        }

        // 5.) Create the Polygon and set needsSort to `true`as we expect it to be convex and
        //     sortable
        return this.createPolygon(coordinates, true);
    }

    private getWorldConers(projection: Projection): TileCorners | undefined {
        if (projection.type !== ProjectionType.Planar) {
            return;
        }
        const worldBox = projection.worldExtent(0, 0);
        return {
            sw: worldBox.min as Vector3,
            se: new Vector3(worldBox.max.x, worldBox.min.y, 0),
            nw: new Vector3(worldBox.min.x, worldBox.max.y, 0),
            ne: worldBox.max as Vector3
        };
    }

    private addNDCRayIntersection(
        ndcPoints: Array<[number, number]>,
        geoPolygon: GeoCoordinates[]
    ) {
        ndcPoints.forEach(corner => {
            const intersection = MapViewUtils.rayCastWorldCoordinates(
                { camera: this.m_camera, projection: this.m_projection },
                corner[0],
                corner[1]
            );
            if (intersection) {
                this.validateAndAddToGeoPolygon(intersection, geoPolygon);
            }
        });
    }

    private addHorizonIntersection(geoPolygon: GeoCoordinates[]) {
        if (this.m_projection.type === ProjectionType.Planar) {
            const verticalHorizonPosition = this.getVerticalHorizonPositionInNDC();
            if (!verticalHorizonPosition) {
                return;
            }
            this.addNDCRayIntersection(
                [
                    [-1, verticalHorizonPosition], //horizon left
                    [1, verticalHorizonPosition] //horizon right
                ],
                geoPolygon
            );
        }
    }

    private addCanvasCornerIntersection(geoPolygon: GeoCoordinates[]) {
        this.addNDCRayIntersection(
            [
                [-1, -1], //lower left
                [1, -1], //lower right
                [1, 1], //upper right
                [-1, 1] //upper left
            ],
            geoPolygon
        );
    }

    private validateAndAddToGeoPolygon(point: Vector3, geoPolygon: GeoCoordinates[]) {
        if (this.isInVisibleMap(point)) {
            geoPolygon.push(this.m_projection.unprojectPoint(point));
        }
    }

    private isInVisibleMap(point: Vector3): boolean {
        if (this.m_projection.type === ProjectionType.Planar) {
            if (point.y < 0 || point.y > EarthConstants.EQUATORIAL_CIRCUMFERENCE) {
                return false;
            }

            if (
                !this.tileWrappingEnabled &&
                (point.x < 0 || point.x > EarthConstants.EQUATORIAL_CIRCUMFERENCE)
            ) {
                return false;
            }
        }
        return true;
    }

    private addPointInFrustum(point: Vector3, frustum: Frustum, geoPolygon: GeoCoordinates[]) {
        if (frustum.containsPoint(point)) {
            const geoPoint = this.m_projection.unprojectPoint(point);
            geoPoint.altitude = 0;
            geoPolygon.push(geoPoint);
        }
    }

    private addFrustumIntersection(
        edge: Line3 | Ray,
        frustum: Frustum,
        geoPolygon: GeoCoordinates[]
    ) {
        frustum.planes.forEach(plane => {
            let intersection: Vector3 | null | undefined = null;
            const target: Vector3 = new Vector3();
            if (edge instanceof Ray && edge.intersectsPlane(plane)) {
                intersection = edge.intersectPlane(plane, target);
            } else if (edge instanceof Line3 && plane.intersectsLine(edge)) {
                intersection = plane.intersectLine(edge, target);
            }

            if (intersection) {
                //uses this check to fix inaccuracies
                if (MapViewUtils.closeToFrustum(intersection, this.m_camera)) {
                    const geoIntersection = this.m_projection.unprojectPoint(intersection);

                    //correct altitude caused by inaccuracies, due to large numbers to 0
                    geoIntersection.altitude = 0;
                    geoPolygon.push(geoIntersection);
                }
            }
        });
    }

    private getVerticalHorizonPositionInNDC(): number | undefined {
        if (this.m_projection.type !== ProjectionType.Planar) {
            return undefined;
        }

        const bottomMidFarPoint = new Vector3(-1, -1, 1)
            .unproject(this.m_camera)
            .add(new Vector3(1, -1, 1).unproject(this.m_camera))
            .multiplyScalar(0.5);
        const topMidFarPoint = new Vector3(-1, 1, 1)
            .unproject(this.m_camera)
            .add(new Vector3(1, 1, 1).unproject(this.m_camera))
            .multiplyScalar(0.5);
        const farPlaneVerticalCenterLine = new Line3(bottomMidFarPoint, topMidFarPoint);

        const verticalHorizonPosition: Vector3 = new Vector3();
        if (
            !this.m_groundPlane.intersectLine(farPlaneVerticalCenterLine, verticalHorizonPosition)
        ) {
            return undefined;
        }
        return verticalHorizonPosition.project(this.m_camera).y;
    }
}
