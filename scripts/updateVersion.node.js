var fs = require('fs');
var path = require('path');

// Constants
var DIR_PATH = './tests';
var DRIVER = 'LatticeConstantHexagonalEnergy__TD_942334626465_000'

// Obtain New Driver Name
var driver = DRIVER;
var newDriver = driver.replace(/\d{3}$/, '001');

// Obtain All Nodes in DIR_PATH
var dirs = fs.readdirSync(DIR_PATH);

for(var i = 0; i < dirs.length; i++) {
    var dir = dirs[i];
    if(/__TE_\d{12}_\d{3}$/.test(dir) == false) {
        // Ignore those that are not tests
        continue;
    }
    
    // Obtain New Test Name
    var test = dir;
    var newTest = test.replace(/\d{3}$/, '001');
    
    // Replace Directory Version
    fs.renameSync(
        path.join(DIR_PATH, test),
        path.join(DIR_PATH, newTest)
    );
    
    var files = fs.readdirSync(path.join(DIR_PATH, newTest));
    for(var j = 0; j < files.length; j++) {
        var file = path.join(DIR_PATH, newTest, files[j]);
        var content = fs.readFileSync(file, 'utf-8');
        var oldContent = content;
        // Replace driver version
        content = content.replace(driver, newDriver);
        // Replace test version
        content = content.replace(test, newTest);
        if(oldContent != content) {
            // Write Back Content
            fs.writeFileSync(file, content);            
            console.log('Updated: ' + file);
            console.log(content);
        }
    }
}