/*
 * Copyright (C) 2017-2020 HERE Europe B.V.
 * Licensed under Apache 2.0, see full license in LICENSE
 * SPDX-License-Identifier: Apache-2.0
 */
import { GeoBox, GeoCoordinates, TileKey, webMercatorTilingScheme } from "@here/harp-geoutils";
import { GeoCoordLike } from "@here/harp-geoutils/lib/coordinates/GeoCoordLike";
import { LongPressHandler, MapControls, MapControlsUI } from "@here/harp-map-controls";
import { CopyrightElementHandler, MapView } from "@here/harp-mapview";
import { DataProvider } from "@here/harp-mapview-decoder";
import { VectorTileDataSource } from "@here/harp-vectortile-datasource";
import { VectorTileDecoder } from "@here/harp-vectortile-datasource/index-worker";
import { Box3, Color } from "three";

import { apikey } from "../config";

export namespace MapObjectsExamples {
    interface Marker {
        coordinates: GeoCoordLike;
        color?: string;
        minZoomLevel?: number;
        maxZoomLevel?: number;
    }

    interface Polygon {
        coordinates: GeoCoordLike[][];
        color?: string;
        minZoomLevel?: number;
        maxZoomLevel?: number;
    }

    function geoBoxFromCoords(coords: GeoCoordLike[]): GeoBox | undefined {
        if (!Array.isArray(coords) || coords.length === 0) {
            return undefined;
        }
        const geoBox = new GeoBox(
            GeoCoordinates.fromObject(coords[0]),
            GeoCoordinates.fromObject(coords[0])
        );
        const reduce = (geoBox: GeoBox, coord: GeoCoordLike) => {
            geoBox.growToContain(GeoCoordinates.fromObject(coord));
            return geoBox;
        };
        return coords.reduce(reduce, geoBox);
    }

    class MapObjectsDataProvider extends DataProvider {
        readonly markers: Marker[] = [];
        readonly polygons: Polygon[] = [];

        ready(): boolean {
            return true;
        }

        async getTile(
            tileKey: TileKey,
            abortSignal?: AbortSignal
        ): Promise<{} | ArrayBuffer | SharedArrayBuffer> {
            const geoBox = webMercatorTilingScheme.getGeoBox(tileKey);

            const tiledMarkers = this.markers.filter(marker =>
                geoBox.contains(GeoCoordinates.fromObject(marker.coordinates))
            );

            const pointFeatures = tiledMarkers.map(marker => {
                const { latitude, longitude } = GeoCoordinates.fromObject(marker.coordinates);
                return {
                    type: "Feature",
                    properties: {
                        lat: latitude.toFixed(2),
                        lng: longitude.toFixed(2),
                        color: marker.color ?? "rgb(255,0,0)",
                        minzoom: marker.minZoomLevel ?? 1,
                        maxzoom: marker.maxZoomLevel ?? 20
                    },
                    geometry: {
                        type: "Point",
                        coordinates: [longitude, latitude]
                    }
                };
            });

            const tileBounds = new Box3();
            webMercatorTilingScheme.getWorldBox(tileKey, tileBounds);
            const polygonBounds = new Box3();
            const interestingPolygons = this.polygons.filter(polygon => {
                const geoBox = geoBoxFromCoords(polygon.coordinates[0]);
                if (geoBox === undefined) {
                    return false;
                }
                webMercatorTilingScheme.projection.projectBox(geoBox, polygonBounds);
                return tileBounds.intersect(polygonBounds);
            });

            const polygonFeatures = interestingPolygons.map(polygon => {
                const coordinates = polygon.coordinates.map(coords =>
                    coords.map(coord => GeoCoordinates.fromObject(coord).toGeoPoint())
                );
                return {
                    type: "Feature",
                    properties: {
                        color: polygon.color ?? "rgb(255,0,0)",
                        minzoom: polygon.minZoomLevel ?? 1,
                        maxzoom: polygon.maxZoomLevel ?? 20
                    },
                    geometry: {
                        type: "Polygon",
                        coordinates
                    }
                };
            });

            return {
                type: "FeatureCollection",
                features: [...pointFeatures, ...polygonFeatures]
            };
        }

        protected async connect() {
            // nothing to do
        }

        protected dispose(): void {
            // nothign to do
        }
    }

    class MapObjects extends VectorTileDataSource {
        private readonly m_dataProvider: MapObjectsDataProvider;

        get markers(): readonly Marker[] {
            return this.m_dataProvider.markers;
        }

        get polygons(): readonly Polygon[] {
            return this.m_dataProvider.polygons;
        }

        constructor() {
            super({
                dataProvider: new MapObjectsDataProvider(),
                styleSetName: "map-objects",
                decoder: new VectorTileDecoder()
            });

            this.m_dataProvider = this.dataProvider() as MapObjectsDataProvider;
        }

        clearMarkers() {
            this.m_dataProvider.markers.length = 0;
            this.mapView.clearTileCache(this.name);
            this.requestUpdate();
        }

        clearPolygons() {
            this.m_dataProvider.polygons.length = 0;
            this.mapView.clearTileCache(this.name);
            this.requestUpdate();
        }

        addMarker(marker: Marker) {
            this.m_dataProvider.markers.push(marker);

            this.mapView.clearTileCache(this.name, tile =>
                tile.geoBox.contains(GeoCoordinates.fromObject(marker.coordinates))
            );

            this.requestUpdate();
        }

        addPolygon(polygon: Polygon) {
            this.m_dataProvider.polygons.push(polygon);

            const geoBox = geoBoxFromCoords(polygon.coordinates[0]);

            if (geoBox === undefined) {
                return;
            }

            const worldBox = new Box3();
            webMercatorTilingScheme.projection.projectBox(geoBox, worldBox);

            const tileBounds = new Box3();
            this.mapView.clearTileCache(this.name, tile => {
                webMercatorTilingScheme.projection.projectBox(tile.geoBox, tileBounds);
                return tileBounds.intersectsBox(worldBox);
            });

            this.requestUpdate();
        }
    }

    async function main(id: string = "mapCanvas") {
        const canvas = document.getElementById(id) as HTMLCanvasElement;

        const NY = new GeoCoordinates(40.707, -74.01);
        const map = new MapView({
            canvas,
            theme: "resources/berlin_tilezen_base.json",
            target: NY,
            zoomLevel: 16.1
        });

        CopyrightElementHandler.install("copyrightNotice", map);

        const mapControls = new MapControls(map);
        mapControls.maxTiltAngle = 50;

        const ui = new MapControlsUI(mapControls, { zoomLevel: "input" });
        canvas.parentElement!.appendChild(ui.domElement);

        map.resize(window.innerWidth, window.innerHeight);

        window.addEventListener("resize", () => {
            map.resize(window.innerWidth, window.innerHeight);
        });

        const omvDataSource = new VectorTileDataSource({
            baseUrl: "https://vector.hereapi.com/v2/vectortiles/base/mc",
            authenticationCode: apikey
        });

        map.addDataSource(omvDataSource);

        const mapObjects = new MapObjects();

        await map.addDataSource(mapObjects);

        mapObjects.setStyleSet([
            {
                styleSet: "map-objects",
                when: ["==", ["geometry-type"], "Polygon"],
                technique: "fill",
                color: ["get", "color"],
                renderOrder: ["number", ["get", "sort-rank"], 100000],
                minZoomLevel: ["get", "minzoom"],
                maxZoomLevel: ["get", "maxzoom"]
            },
            {
                styleSet: "map-objects",
                when: ["==", ["geometry-type"], "Point"],
                technique: "circles",
                size: 50,
                color: ["get", "color"],
                renderOrder: ["number", ["get", "sort-rank"], 100000],
                minZoomLevel: ["get", "minzoom"],
                maxZoomLevel: ["get", "maxzoom"]
            }
        ]);

        const handler = new LongPressHandler(canvas, event => {
            const coordinates = map.getGeoCoordinatesAt(event.pageX, event.pageY);

            if (!coordinates) {
                return;
            }

            const randomColor = new Color(Math.random(), Math.random(), Math.random());

            const color = `#${randomColor.getHexString()}`;

            mapObjects.addMarker({
                coordinates,
                color,
                minZoomLevel: 13,
                maxZoomLevel: 16
            });

            if (mapObjects.markers.length % 3 === 0) {
                const N = mapObjects.markers.length / 3;
                mapObjects.clearPolygons();
                for (let i = 0; i < N; ++i) {
                    mapObjects.addPolygon({
                        coordinates: [
                            [
                                mapObjects.markers[i * 3].coordinates,
                                mapObjects.markers[i * 3 + 1].coordinates,
                                mapObjects.markers[i * 3 + 2].coordinates
                            ]
                        ],
                        color: "rgba(255,255,0,0.5)",
                        minZoomLevel: 13,
                        maxZoomLevel: 16
                    });
                }
                return;
            }
        });

        handler.timeout = 200;
    }

    // eslint-disable-next-line no-console
    main().catch(console.error);
}
