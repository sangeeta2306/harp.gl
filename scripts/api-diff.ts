import { ApiDeclaredItem, ApiItem, ApiModel } from "@microsoft/api-extractor-model";
import * as fs from "fs";
import * as glob from "glob";
import { EOL } from "os";
import * as path from "path";

function traverse(apiItem: ApiItem, depth: number, stream: fs.WriteStream) {
    const indent = " ".repeat(depth * 2);
    const kind = apiItem.kind;
    if (apiItem instanceof ApiDeclaredItem) {
        const text = apiItem.excerpt.text.split("\n").join(`${EOL}${indent}`);
        stream.write(`${indent}${text}${EOL}`);
    } else {
        stream.write(`${indent}${apiItem.displayName}: ${kind}${EOL}`);
    }
    apiItem.members.forEach(childItem => {
        traverse(childItem, depth + 1, stream);
    });
    if (apiItem.members.length > 0) {
        stream.write(EOL);
    }
}

async function main() {
    const stream = fs.createWriteStream("API.txt");
    const moduleAPIs = glob.sync("input/*.json");
    moduleAPIs.forEach(fileName => {
        // const outFileName = path.basename(fileName).replace(".json", ".txt");
        // console.log(`Writing ${outFileName}`);
        const apiModel = new ApiModel();
        apiModel.loadPackage(fileName);
        traverse(apiModel, 0, stream);
    });
    stream.close();
}

main().catch(console.error);
