#ifndef XMLSERIALIZER_H
#define XMLSERIALIZER_H

#include <QXmlStreamWriter>

#include "graphserialization_global.h"
#include "serializer.h"

class GRAPHSERIALIZATION_EXPORT XMLSerializer : public Serializer
{
public:
    XMLSerializer();

    void save(const Serializable &serializable, const QString &filepath, const QString &rootName = "") override;

    void load(Serializable &serializable, const QString &filepath) override;

private:

    void writeVariantToStream(const QString& nodeName, const QVariant& variant, QXmlStreamWriter& stream);
    void writeVariantValueToStream(const QVariant& variant, QXmlStreamWriter& stream);
    void writeVariantListToStream(const QVariant& variantList, QXmlStreamWriter& stream);
    void writeVariantMapToStream(const QVariant& variantMap, QXmlStreamWriter& stream);


    QVariant readVariantFromStream(QXmlStreamReader& stream);
    QVariant readVariantValueFromStream(QXmlStreamReader& stream);
    QVariant readVariantListFromStream(QXmlStreamReader& stream);
    QVariant readVariantMapFromStream(QXmlStreamReader& stream);
};


#endif // XMLSERIALIZER_H
