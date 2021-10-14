module Units

export second, minute, hour, meter, kilometer,
       seconds, minutes, hours, meters, kilometers,
       KiB, MiB, GiB, TiB

const second = 1.0
const seconds = second
const minute = 60seconds
const minutes = minute
const hour = 60minutes
const hours = hour

const meter = 1.0
const meters = meter
const kilometer = 1000meters
const kilometers = kilometer

const KiB = 1024.0
const MiB = 1024KiB
const GiB = 1024MiB
const TiB = 1024GiB

end